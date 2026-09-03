# deconstructSigs: panel opportunity correction AND SATS's de novo signature filter
#
# The two corrections this project cares about, applied together:
#
#   1. the panel correction -- each channel is scaled by the assay's own target
#      territory instead of assuming the exome (deconstructsigs_panel_assignment.R)
#   2. the de novo filter -- each tumour is fitted only against the signatures
#      SATS's discovery step selected for its cancer cohort
#
# They compose without interacting. deconstructSigs takes the opportunity
# correction into the DATA (a 32-value ratio applied to the counts, by assay),
# while the filter restricts the SIGNATURES (by cancer type). Neither touches the
# other, so no per-(cohort x assay) fitting is needed here -- unlike sigminer and
# SigProfiler, which take the correction into the signatures and therefore need
# both keys at once.
#
# There is also a consistency argument for combining them: SATS ran its de novo
# WITH the panel L, so the signature lists were chosen under a panel opportunity
# model. Fitting them against exome-corrected data would take the candidate set
# from one model and apply it under another.
#
# No de novo is repeated here; the lists come from the cohort_signatures.tsv the
# SATS by-cohort script writes.
#
# THE ONE TRANSLATION
#
# SATS maps to RefTMB, where APOBEC is the merged column SBS2_13; the
# deconstructSigs reference carries SBS2 and SBS13 separately, so SBS2_13 expands
# to both. The other 84 names transfer verbatim.
#
# Output: <results_dir>/<out_name>/deconstructSigs_assignment.csv
#         one row per tumour, one column per COSMIC v3.4 signature, zero where
#         that signature was not in the tumour's cohort list.

library(readr)
library(deconstructSigs)
library(dplyr)
library(tibble)
library(data.table)
library(rlang)
setwd("C:/Users/Josep/Desktop/AACR_local/consensus-sigs-main/consensus-sigs-main/scripts")
source("C:/Users/Josep/Desktop/AACR_local/consensus-sigs-main/consensus-sigs-main/scripts/utils/panel_opportunity_utils.R")

## helpers 

expand_apobec <- function(x) {
  if (!"SBS2_13" %in% x) return(x)
  unique(c(setdiff(x, "SBS2_13"), "SBS2", "SBS13"))
}

read_cohort_signatures <- function(path, valid) {
  if (!file.exists(path))
    stop("cannot find '", path, "'. Run the SATS by-cohort script first; it ",
         "writes cohort_signatures.tsv beside its assignment file.")
  d <- data.table::fread(path)
  if (!all(c("cohort", "signature") %in% names(d)))
    stop("cohort_signatures.tsv must have columns 'cohort' and 'signature'")
  out <- lapply(split(as.character(d$signature), as.character(d$cohort)),
                function(s) intersect(expand_apobec(unique(s)), valid))
  out[lengths(out) > 0]
}

## driver 

run_deconstructsigs_panel_filtered <- function(
    maf, results_dir, genomic_info, clinical,
    signatures_from,
    cohort_column = "CANCER_TYPE",
    unmatched     = c("skip", "all"),
    cores         = 1L,
    out_name      = "deconstructSigs_PanelFiltered_Cosmic_V3.4") {

  unmatched <- match.arg(unmatched)
  out_dir <- file.path(results_dir, out_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  maf_table <- readr::read_tsv(maf, show_col_types = FALSE)
  cosmic <- readRDS("./resources/deconstruct_sigs/signatures_SBS_cosmic.v3.4.rds")
  all_sigs <- rownames(cosmic)

  ## ---- correction 1: the panel ratio, into the data ----------------------
  message("building panel opportunity from ", genomic_info)
  intervals  <- read_panel_intervals(genomic_info)
  sample_ids <- unique(maf_table$Tumor_Sample_Barcode)
  groups     <- panel_groups(sample_ids, intervals, clinical)
  message("  ", length(groups), " assay group(s)")

  ratio_of <- lapply(names(groups), function(g)
    panel2genome_ratio(panel_trinuc_counts(panel_intervals_for(intervals, g))))
  names(ratio_of) <- names(groups)

  # chr must be a FACTOR: mut.to.sigs.input() renames "5" to "chr5" with
  # levels(mut[, chr]) <- sub(...), which only rewrites values on a factor.
  all_mutations <- maf_table %>%
    dplyr::select(Sample = Tumor_Sample_Barcode, chr = Chromosome,
                  pos = Start_Position, ref = Reference_Allele,
                  alt = Tumor_Seq_Allele2) %>%
    dplyr::mutate(chr = as.factor(.data$chr))

  sigs_input_all <- mut.to.sigs.input(
    mut.ref = as.data.frame(all_mutations), sample.id = "Sample",
    chr = "chr", pos = "pos", ref = "ref", alt = "alt")

  for (g in names(groups)) {
    rows <- intersect(groups[[g]], rownames(sigs_input_all))
    if (length(rows))
      sigs_input_all[rows, ] <- apply_panel_ratio(
        sigs_input_all[rows, , drop = FALSE], ratio_of[[g]])
  }
  message("  counts corrected: ", nrow(sigs_input_all), " tumours x ",
          ncol(sigs_input_all), " channels")

  ## ---- correction 2: the de novo filter, into the signatures -------------
  cohort_sigs <- read_cohort_signatures(signatures_from, all_sigs)
  message("signature lists read for ", length(cohort_sigs), " cohorts")

  cl <- data.table::fread(clinical, sep = "\t", quote = "", skip = "PATIENT_ID",
                          showProgress = FALSE)
  if (!all(c("SAMPLE_ID", cohort_column) %in% names(cl)))
    stop("clinical table needs SAMPLE_ID and ", cohort_column)
  ct <- setNames(as.character(cl[[cohort_column]]), cl$SAMPLE_ID)[rownames(sigs_input_all)]
  has_list <- !is.na(ct) & ct %in% names(cohort_sigs)
  message("  tumours in a cohort with a list: ", sum(has_list), " of ",
          nrow(sigs_input_all))
  if (any(!has_list))
    message("  ", sum(!has_list), " without a list (",
            if (unmatched == "skip") "skipped" else "all 86 signatures", ")")

  jobs <- split(rownames(sigs_input_all)[has_list], ct[has_list])
  if (unmatched == "all" && any(!has_list))
    jobs[["__all_signatures__"]] <- rownames(sigs_input_all)[!has_list]

  ## ---- fit ----------------------------------------------------------------
  # tri.counts.method = "default": the panel correction is already in the counts,
  # so any further reweighting here would correct twice.
  W_all <- matrix(0, nrow(sigs_input_all), length(all_sigs),
                  dimnames = list(rownames(sigs_input_all), all_sigs))
  summ <- list()
  for (cname in names(jobs)) {
    rows <- jobs[[cname]]
    sigs <- if (identical(cname, "__all_signatures__")) all_sigs else cohort_sigs[[cname]]
    message("\n--- ", cname, ": ", length(rows), " tumours, ", length(sigs),
            " signatures: ", paste(sigs, collapse = ", "))
    ref_sub <- cosmic[sigs, , drop = FALSE]
    one <- function(s) tryCatch(
      whichSignatures(tumor.ref = sigs_input_all[s, , drop = FALSE],
                      signatures.ref = ref_sub,
                      tri.counts.method = "default",
                      contexts.needed = TRUE)$weights,
      error = function(e) NULL)
    got <- if (cores > 1L) {
      clus <- parallel::makeCluster(cores)
      on.exit(parallel::stopCluster(clus), add = TRUE)
      parallel::clusterEvalQ(clus, library(deconstructSigs))
      parallel::clusterExport(clus, c("sigs_input_all", "ref_sub"), envir = environment())
      parallel::parLapply(clus, rows, one)
    } else lapply(rows, one)
    names(got) <- rows
    ok <- !vapply(got, is.null, logical(1))
    for (s in rows[ok]) W_all[s, colnames(got[[s]])] <- as.numeric(got[[s]][1, ])
    if (any(!ok)) message("  ", sum(!ok), " tumour(s) produced no fit")
    summ[[cname]] <- data.table(cohort = cname, n_tumours = length(rows),
                                n_fitted = sum(ok), n_signatures = length(sigs),
                                signatures = paste(sigs, collapse = ","))
  }

  keep <- if (unmatched == "skip") rownames(W_all)[has_list] else rownames(W_all)
  out <- tibble::rownames_to_column(as.data.frame(W_all[keep, , drop = FALSE]),
                                    var = "sample_name")
  readr::write_csv(out, file.path(out_dir, "deconstructSigs_assignment.csv"))
  S <- data.table::rbindlist(summ, fill = TRUE)
  data.table::fwrite(S, file.path(out_dir, "cohort_summary.tsv"), sep = "\t")

  cat("\n================ SUMMARY ================\n")
  print(S[, .(cohort, n_tumours, n_fitted, n_signatures)], row.names = FALSE)
  cat("\ntumours written: ", nrow(out),
      "\nmean signatures per cohort: ", round(mean(S$n_signatures), 1),
      " (versus 86 unfiltered)",
      "\nwritten: ", file.path(out_dir, "deconstructSigs_assignment.csv"), "\n", sep = "")
  invisible(list(out = out, summary = S))
}


run_deconstructsigs_panel_filtered(
  maf             = "../../../data/input/19.0_public/data_mutations_extended.txt",
  results_dir     = "test_sats/",
  genomic_info    = "../../../data/input/19.0_public/genomic_information.txt",
  clinical        = "../../../data/input/19.0_public/data_clinical_sample.txt",
  signatures_from = "test_sats/SATS_ByCohort_Cosmic_V3.4/cohort_signatures.tsv",
  cohort_column   = "CANCER_TYPE",
  cores           = 8L)
