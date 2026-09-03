# sigminer + SigProfiler: panel opportunity correction AND SATS's de novo filter
#
# Both corrections at once:
#
#   1. the panel correction -- COSMIC rescaled to the assay's own target
#      territory, W_panel proportional to W_genome * L_panel / L_genome
#   2. the de novo filter -- only the signatures SATS's discovery selected for
#      that tumour's cancer cohort
#
# WHY THIS ONE NEEDS A CELL LOOP AND deconstructSigs DOES NOT
#
# Neither tool has an opportunity term, so for both of them the panel correction
# goes into the SIGNATURES. The de novo filter also acts on the signatures. So
# the signature matrix handed to each fit must be rescaled for the assay AND
# subset to the cohort at the same time -- which means fitting per
# (cancer type x assay) CELL, not per cohort.
#
# In GENIE 19.0 that is 3,421 cells that actually occur (of 18,704 possible),
# 1,434 with >= 10 tumours covering 97.8% of tumours. deconstructSigs escapes
# this because it takes the correction into the DATA, so the two keys never meet.
#
# The saving grace: rescaling and subsetting COMMUTE.
# rescale_signatures_to_panel() normalises each column by its own sum, so the
# columns are independent. Verified on MSK-IMPACT505: rescale-then-subset versus
# subset-then-rescale gives max absolute difference 0. So the expensive part --
# building each assay's trinucleotide territory and rescaling all 86 signatures
# -- is done ONCE PER ASSAY and cached; each cell then just takes columns.
#
# Cells are checkpointed individually and a re-run skips completed ones, because
# a full pass is thousands of fits.
#
# THE ONE TRANSLATION
#
# SATS maps to RefTMB, where APOBEC is the merged column SBS2_13; both references
# here carry SBS2 and SBS13 separately, so SBS2_13 expands to both.
#
# A NOTE ON SCALE
#
# Unlike the unfiltered arms, both tools here are on the SAME scale: the
# panel-rescaled catalogue. The asymmetry in the unfiltered arms (sigminer
# genome, SigProfiler exome) exists because sig_fit() silently discards
# `exome = TRUE` while cosmic_fit() honours it. Supplying an explicit rescaled
# matrix to both removes that discrepancy -- which is the point of the panel
# correction, not a side effect to apologise for.
#
# Output:
#   <results_dir>/<out_sigminer>/SigMiner_assignment.csv
#   <results_dir>/<out_sigprofiler>/<cell>/... plus combined_Activities.txt

library(dplyr)
library(maftools)
library(reticulate)
library(SigProfilerAssignmentR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(quadprog)
library(sigminer)
library(data.table)
setwd("C:/Users/Josep/Desktop/AACR_local/consensus-sigs-main/consensus-sigs-main/scripts")
source("./utils/panel_opportunity_utils.R")

python_path <- Sys.getenv("RETICULATE_PYTHON", unset = NA)
if (!is.na(python_path) && nzchar(python_path)) use_python(python_path, required = TRUE)

##helpers

expand_apobec <- function(x) {
  if (!"SBS2_13" %in% x) return(x)
  unique(c(setdiff(x, "SBS2_13"), "SBS2", "SBS13"))
}

read_cohort_signatures <- function(path, valid) {
  if (!file.exists(path))
    stop("cannot find '", path, "'. Run the SATS by-cohort script first.")
  d <- data.table::fread(path)
  if (!all(c("cohort", "signature") %in% names(d)))
    stop("cohort_signatures.tsv must have columns 'cohort' and 'signature'")
  out <- lapply(split(as.character(d$signature), as.character(d$cohort)),
                function(s) intersect(expand_apobec(unique(s)), valid))
  out[lengths(out) > 0]
}

# SigProfiler reads the channel order from the file it is given. Write the
# signature matrix and the sample matrix in the SAME order so the two cannot
# disagree, and assert the round trip.
write_spa_pair <- function(W, counts, channels, sig_file, mat_file) {
  ord <- sort(channels)
  sig <- data.frame(MutationType = ord, W[ord, , drop = FALSE], check.names = FALSE)
  data.table::fwrite(sig, sig_file, sep = "\t", quote = FALSE, row.names = FALSE)
  back <- data.table::fread(sig_file, sep = "\t", header = TRUE)
  if (nrow(back) != 96L || !identical(as.character(back$MutationType), ord))
    stop("the custom signature file lost its 96-channel ascending order")
  mat <- data.frame(MutationType = ord, counts[ord, , drop = FALSE], check.names = FALSE)
  data.table::fwrite(mat, mat_file, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(TRUE)
}

##driver

run_sigminer_sigprofiler_panel_filtered <- function(
    maf, results_dir, genomic_info, clinical,
    signatures_from,
    cohort_column    = "CANCER_TYPE",
    min_cell_samples = 1L,
    run_sigprofiler  = TRUE,
    resume           = TRUE,
    out_sigminer     = "SigMiner_PanelFiltered_Cosmic_V3.4",
    out_sigprofiler  = "SigProfiler_PanelFiltered_Cosmic_V3.4",
    tmp_dir          = "./tmp/panel_filtered") {

  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  sm_dir <- file.path(results_dir, out_sigminer)
  sp_dir <- file.path(results_dir, out_sigprofiler)
  ckpt   <- file.path(results_dir, out_sigminer, "checkpoints")
  for (d in c(sm_dir, sp_dir, ckpt)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  ## ---- counts, once -------------------------------------------------------
  maf_obj <- sigminer::read_maf(maf = maf)
  tally <- sigminer::sig_tally(maf_obj, ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                               genome_build = "hg19", useSyn = TRUE)
  catalogue <- t(tally$nmf_matrix)                       # 96 x N
  channels  <- rownames(catalogue)
  message("counts: ", nrow(catalogue), " channels x ", ncol(catalogue), " tumours")

  cosmic <- readRDS("./resources/deconstruct_sigs/signatures_SBS_cosmic.v3.4.rds")
  if (!setequal(colnames(cosmic), channels))
    stop("the tally and the COSMIC reference do not share the same 96 channels")
  all_sigs <- rownames(cosmic)

  ## ---- keys: assay (for the correction) and cohort (for the filter) -------
  intervals <- read_panel_intervals(genomic_info)
  assay_groups <- panel_groups(colnames(catalogue), intervals, clinical)
  assay_of <- setNames(rep(names(assay_groups), lengths(assay_groups)),
                       unlist(assay_groups))

  cohort_sigs <- read_cohort_signatures(signatures_from, all_sigs)
  cl <- data.table::fread(clinical, sep = "\t", quote = "", skip = "PATIENT_ID",
                          showProgress = FALSE)
  if (!all(c("SAMPLE_ID", cohort_column) %in% names(cl)))
    stop("clinical table needs SAMPLE_ID and ", cohort_column)
  ct <- setNames(as.character(cl[[cohort_column]]), cl$SAMPLE_ID)[colnames(catalogue)]

  usable <- !is.na(ct) & ct %in% names(cohort_sigs) & !is.na(assay_of[colnames(catalogue)])
  message("signature lists: ", length(cohort_sigs), " cohorts; assays: ",
          length(assay_groups))
  message("  tumours with both a cohort list and an assay: ", sum(usable),
          " of ", ncol(catalogue))

  # Cells are held as a table of (cohort, assay, sample), never as a pasted
  # string key. Cohort names contain spaces and commas, assay names contain
  # hyphens and dots, so any separator risks colliding with the data and
  # splitting a composite key back apart is a silent-corruption bug waiting to
  # happen.
  cell_dt <- data.table::data.table(
    sample = colnames(catalogue)[usable],
    cohort = unname(ct[usable]),
    assay  = unname(assay_of[colnames(catalogue)[usable]]))
  cells <- cell_dt[, .(samples = list(sample), n = .N), by = .(cohort, assay)]
  cells <- cells[n >= min_cell_samples][order(-n)]
  message("  (cohort x assay) cells to fit: ", nrow(cells),
          " (min ", min_cell_samples, " tumours)")

  ## ---- rescale ONCE PER ASSAY --------------------------------------------
  # This is the expensive half: each assay's target territory and a rescale of
  # all 86 signatures. Because rescaling and subsetting commute, each cell then
  # only has to take columns from the cached matrix.
  needed_assays <- unique(cells$assay)
  message("  rescaling COSMIC for ", length(needed_assays), " assay(s)")
  W_of <- new.env(parent = emptyenv())
  for (a in needed_assays) {
    tri <- panel_trinuc_counts(panel_intervals_for(intervals, a))
    assign(a, rescale_signatures_to_panel(cosmic, tri, channels = channels), envir = W_of)
  }

  ## ---- per cell ------------------------------------------------------------
  sm_parts <- list(); summ <- list()
  for (i in seq_len(nrow(cells))) {
    cname <- cells$cohort[i]; aname <- cells$assay[i]
    cols  <- cells$samples[[i]]
    key   <- paste(cname, aname, sep = " | ")
    slug <- gsub("[^A-Za-z0-9]+", "_", paste(cname, aname, sep = "__"))
    cfile <- file.path(ckpt, paste0(slug, ".rds"))

    if (resume && file.exists(cfile)) {
      got <- readRDS(cfile)
      summ[[key]] <- got$rec
      if (!is.null(got$sm)) sm_parts[[key]] <- got$sm
      next
    }

    sigs <- intersect(cohort_sigs[[cname]], all_sigs)
    W <- get(aname, envir = W_of)[, sigs, drop = FALSE]   # rescaled AND filtered
    rec <- data.table(cohort = cname, assay = aname, n_tumours = length(cols),
                      n_signatures = length(sigs), sigminer = FALSE,
                      sigprofiler = FALSE)
    message("\n--- ", cname, " | ", aname, ": ", length(cols), " tumours, ",
            length(sigs), " signatures")

    part <- tryCatch(
      sigminer::sig_fit(catalogue_matrix = catalogue[, cols, drop = FALSE],
                        sig = W, return_class = "data.table"),
      error = function(e) { message("  sigminer FAILED: ", conditionMessage(e)); NULL })
    if (!is.null(part)) { sm_parts[[key]] <- part; rec$sigminer <- TRUE }

    if (run_sigprofiler) {
      sig_file <- file.path(tmp_dir, paste0("sig_", slug, ".txt"))
      mat_file <- file.path(tmp_dir, paste0("mat_", slug, ".txt"))
      rec$sigprofiler <- tryCatch({
        write_spa_pair(W, catalogue[, cols, drop = FALSE], channels, sig_file, mat_file)
        cosmic_fit(mat_file, output = file.path(sp_dir, slug),
          input_type = "matrix", context_type = "96", collapse_to_SBS96 = TRUE,
          cosmic_version = 3.4,
          exome = FALSE,               # the signatures already carry the correction
          genome_build = "GRCh37",
          signature_database = sig_file, exclude_signature_subgroups = NULL,
          export_probabilities = TRUE, export_probabilities_per_mutation = FALSE,
          make_plots = FALSE, sample_reconstruction_plots = FALSE, verbose = FALSE,
          nnls_add_penalty = 0.05, nnls_remove_penalty = 0.0,
          initial_remove_penalty = 0.0)
        TRUE
      }, error = function(e) { message("  SigProfiler FAILED: ", conditionMessage(e)); FALSE })
    }
    summ[[key]] <- rec
    saveRDS(list(rec = rec, sm = sm_parts[[key]]), cfile)
  }

  ## ---- combine -------------------------------------------------------------
  if (length(sm_parts)) {
    sm <- data.table::rbindlist(sm_parts, use.names = TRUE, fill = TRUE)
    for (j in names(sm)) if (is.numeric(sm[[j]]))
      data.table::set(sm, which(is.na(sm[[j]])), j, 0)
    miss <- setdiff(all_sigs, names(sm))
    if (length(miss)) sm[, (miss) := 0]
    data.table::setcolorder(sm, c(setdiff(names(sm), all_sigs), all_sigs))
    data.table::fwrite(sm, file.path(sm_dir, "SigMiner_assignment.csv"), sep = ",")
    message("\nsigminer: ", nrow(sm), " tumours x ", length(all_sigs), " signatures")
  } else message("\nsigminer: nothing fitted")

  acts <- list.files(sp_dir, pattern = "Assignment_Solution_Activities.txt$",
                     recursive = TRUE, full.names = TRUE)
  if (length(acts)) {
    sp <- data.table::rbindlist(lapply(acts, data.table::fread),
                                use.names = TRUE, fill = TRUE)
    for (j in names(sp)) if (is.numeric(sp[[j]]))
      data.table::set(sp, which(is.na(sp[[j]])), j, 0)
    miss <- setdiff(all_sigs, names(sp))
    if (length(miss)) sp[, (miss) := 0]
    data.table::fwrite(sp, file.path(sp_dir, "combined_Activities.txt"), sep = "\t")
    message("SigProfiler: ", nrow(sp), " tumours from ", length(acts), " cell runs")
  } else message("SigProfiler: no activity files produced")

  S <- data.table::rbindlist(summ, fill = TRUE)
  data.table::fwrite(S, file.path(results_dir, "panel_filtered_cell_summary.tsv"),
                     sep = "\t")
  cat("\n================ SUMMARY ================\n")
  cat("cells fitted: ", nrow(S), "\n",
      "  sigminer ok: ", sum(S$sigminer), "\n",
      "  SigProfiler ok: ", sum(S$sigprofiler), "\n",
      "mean signatures per cell: ", round(mean(S$n_signatures), 1),
      " (versus 86 unfiltered)\n", sep = "")
  invisible(list(summary = S))
}


run_sigminer_sigprofiler_panel_filtered(
  maf             = "../../../data/input/19.0_public/data_mutations_extended.txt",
  results_dir     = "test_sats/",
  genomic_info    = "../../../data/input/19.0_public/genomic_information.txt",
  clinical        = "../../../data/input/19.0_public/data_clinical_sample.txt",
  signatures_from = "test_sats/SATS_ByCohort_Cosmic_V3.4/cohort_signatures.tsv",
  cohort_column   = "CANCER_TYPE")
