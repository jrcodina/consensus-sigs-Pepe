# It follows the workflow in the SATS repository
# (github.com/binzhulab/SATS, User_Guide_SATS_v1.0.10.md):
#
#   1. matched mutation-count (V) and panel-context (L) matrices
#   2. de novo signature extraction with signeR, using L as the opportunity
#   3. MappingSignature() onto the TMB-normalised COSMIC catalogue
#   4. EstimateSigActivity() on the mapped subset, then
#      CalculateSignatureBurdens()
#
# Notes on choices that are not obvious:
#
# Output, under <results_dir>/<out_name>/:
#   SATS_assignment_by_cohort.csv   one row per tumour, one column per COSMIC
#                                   v3.4 SBS signature, zero where a cohort did
#                                   not select that signature
#   cohort_summary.tsv              per cohort: n, pools, K, signatures, status
#   cohort_signatures.tsv           long form, cohort x signature
#   checkpoints/<cohort>.rds        per-cohort results, for resume
#

library(data.table)
library(SATS)
library(GenomicRanges)
library(IRanges)
library(Biostrings)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg19)

## ---------------------------------------------------------------- channels ---

# the 96 SBS channels in COSMIC order, e.g. "A[C>A]A"
sbs96_order <- function() {
  subs  <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  flank <- c("A", "C", "G", "T")
  unlist(lapply(subs, function(s)
    as.vector(t(outer(flank, flank, function(a, b) paste0(a, "[", s, "]", b))))))
}

# signeR writes "ACA>A" for what COSMIC writes "A[C>A]A"
signer_to_cosmic <- function(x) {
  if (!all(grepl("^[ACGT]{3}>[ACGT]$", x))) return(NULL)
  tri <- substr(x, 1, 3); alt <- substr(x, 5, 5)
  paste0(substr(tri, 1, 1), "[", substr(tri, 2, 2), ">", alt, "]", substr(tri, 3, 3))
}

## ------------- V and L ---------------------------

# 96 x N mutation count matrix from a MAF-like table.
#
# SATS >= 1.0.10 on GitHub exports GenerateVMatrix() for this; the CRAN build of
# the same version does not (9 exports vs 5). Use it when present, otherwise
# build the matrix here. Both produce counts on the 96 COSMIC channels folded to
# the pyrimidine strand.

build_V <- function(maf_dt, channels, genome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
                    samples = NULL) {
  if ("GenerateVMatrix" %in% getNamespaceExports("SATS")) {
    message("  V: SATS::GenerateVMatrix()")
    V <- as.matrix(SATS::GenerateVMatrix(as.data.frame(maf_dt), Class = "SBS",
                                         ref.genome = "hg19"))
    V <- V[channels, , drop = FALSE]
  } else {
    message("  V: built locally (SATS::GenerateVMatrix not in this build)")
    mut <- maf_dt[Variant_Type == "SNP" &
                  Reference_Allele %in% c("A", "C", "G", "T") &
                  Tumor_Seq_Allele2 %in% c("A", "C", "G", "T") &
                  Reference_Allele != Tumor_Seq_Allele2 &
                  as.character(Chromosome) %in% c(as.character(1:22), "X")]
    if (!nrow(mut)) stop("no usable SNVs")
    gr <- GRanges(paste0("chr", as.character(mut$Chromosome)),
                  IRanges(mut$Start_Position - 1L, mut$Start_Position + 1L))
    sl <- GenomeInfoDb::seqlengths(genome)[as.character(seqnames(gr))]
    inb <- start(gr) >= 1L & end(gr) <= sl
    gr <- gr[inb]; mut <- mut[inb]
    ctx <- as.character(Biostrings::getSeq(genome, gr))
    agree <- substr(ctx, 2, 2) == mut$Reference_Allele
    if (mean(agree) < 0.99)
      warning(sprintf("only %.1f%% of SNVs match hg19 at the reference base",
                      100 * mean(agree)))
    ctx <- ctx[agree]; mut <- mut[agree]
    alt <- mut$Tumor_Seq_Allele2
    pur <- substr(ctx, 2, 2) %in% c("A", "G")
    ctx[pur] <- as.character(Biostrings::reverseComplement(
                  Biostrings::DNAStringSet(ctx[pur])))
    alt[pur] <- chartr("ACGT", "TGCA", alt[pur])
    chan <- paste0(substr(ctx, 1, 1), "[", substr(ctx, 2, 2), ">", alt, "]",
                   substr(ctx, 3, 3))
    mut[, .channel := chan]
    mut <- mut[.channel %in% channels]
    tab <- mut[, .N, by = .(.channel, Tumor_Sample_Barcode)]
    if (is.null(samples)) samples <- sort(unique(mut$Tumor_Sample_Barcode))
    samples <- sort(unique(as.character(samples)))
    V <- matrix(0L, length(channels), length(samples),
                dimnames = list(channels, samples))
    V[cbind(match(tab$.channel, channels),
            match(tab$Tumor_Sample_Barcode, samples))] <- tab$N
  }
  V
}

# 96 x N panel-context matrix, via the SATS API:
#
#   GeneratePanelSize() -> per-assay context counts
#   GenerateLMatrix()   -> per-sample L
#   GenerateLMatrix() requires a data frame with columns PATIENT_ID and
# SEQ_ASSAY_ID; the tumour barcode goes in PATIENT_ID so L is keyed to samples.

build_L <- function(genomic_info, clinical, samples, channels, cache = NULL) {
  if (!is.null(cache) && file.exists(cache)) {
    message("  L: reusing panel context from ", cache)
    Panel_context <- readRDS(cache)
  } else {
    gi <- data.table::fread(genomic_info, showProgress = FALSE)
    need <- c("Chromosome", "Start_Position", "End_Position", "SEQ_ASSAY_ID")
    if (!all(need %in% names(gi)))
      stop("genomic_information must contain: ", paste(need, collapse = ", "))
    gi <- gi[, ..need]
    gi[, Chromosome := sub("^chr", "", as.character(Chromosome))]
    gi <- gi[Chromosome %in% c(as.character(1:22), "X", "Y")]

    # GeneratePanelSize() widens each interval by a base either side to read the
    # trinucleotide context, and does NOT clamp to the chromosome. On GENIE's
    # real interval file that runs off the end and getSeq() aborts with
    #   "trying to load regions beyond the boundaries of non-circular sequence"
    # which SATS then re-reports as "Make sure package BSgenome... is loaded".
    # Clamp to [2, length-1] so the widened interval always stays in bounds.
    genome <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
    sl <- GenomeInfoDb::seqlengths(genome)[paste0("chr", gi$Chromosome)]
    n0 <- nrow(gi)
    gi <- gi[!is.na(sl)]; sl <- sl[!is.na(sl)]
    gi[, Start_Position := pmax(as.integer(Start_Position), 2L)]
    gi[, End_Position   := pmin(as.integer(End_Position), as.integer(sl) - 1L)]
    bad <- gi$End_Position < gi$Start_Position
    if (any(bad))
      message("  L: dropped ", sum(bad), " interval(s) lying outside hg19; ",
              "clamped the rest to chromosome bounds")
    gi <- gi[!bad]
    if (nrow(gi) < n0 - sum(bad))
      message("  L: dropped ", n0 - nrow(gi) - sum(bad),
              " interval(s) on contigs absent from hg19")

    message("  L: SATS::GeneratePanelSize() over ", nrow(gi), " intervals, ",
            uniqueN(gi$SEQ_ASSAY_ID), " assays -- this is the slow step")
    Panel_context <- SATS::GeneratePanelSize(
      genomic_information = as.data.frame(gi), Class = "SBS",
      SBS_order = "COSMIC", ref.genome = "hg19")
    if (!is.null(cache)) saveRDS(Panel_context, cache)
  }

  cl <- data.table::fread(clinical, sep = "\t", quote = "", skip = "PATIENT_ID",
                          showProgress = FALSE)
  if (!all(c("SAMPLE_ID", "SEQ_ASSAY_ID") %in% names(cl)))
    stop("clinical table must contain SAMPLE_ID and SEQ_ASSAY_ID")
  pi <- unique(cl[SAMPLE_ID %in% samples, .(PATIENT_ID = SAMPLE_ID, SEQ_ASSAY_ID)])
  pi <- pi[SEQ_ASSAY_ID %in% colnames(Panel_context)]
  if (!nrow(pi)) stop("no sample maps to an assay present in the panel context")

  L <- as.matrix(SATS::GenerateLMatrix(Panel_context, as.data.frame(pi)))
  L[channels, , drop = FALSE]
}

## ----------------- de novo ------------------------

# signeR de novo extraction, returning 96 x K profiles on COSMIC channel labels.
de_novo <- function(V, L, ranks = 1:12, seed = 1L,
                    testing_burn = 500L, testing_eval = 500L,
                    main_burn = 3000L, main_eval = 1000L,
                    parallelization = "multisession") {

    suppressPackageStartupMessages(library(signeR))

  # signeR's model selection divides by the opportunity; a structural zero there
  # propagates as a non-finite likelihood. Keep the zero in V, floor it in L.
  Lz <- L
  zero <- Lz <= 0
  if (any(zero)) {
    Lz[zero] <- min(Lz[!zero]) * 1e-8
    message("  floored ", sum(zero), " zero opportunity entries for signeR")
  }

  set.seed(seed)
  sr <- signeR(M = t(V), Opport = t(Lz),
               nlim = c(min(ranks), max(ranks)),
               testing_burn = testing_burn, testing_eval = testing_eval,
               main_burn = main_burn, main_eval = main_eval,
               parallelization = parallelization)

  W <- if (!is.null(sr$Phat)) as.matrix(sr$Phat)
       else t(as.matrix(signeR::Median_sign(sr$SignExposures)))
  if (nrow(W) != nrow(V) && ncol(W) == nrow(V)) W <- t(W)
  if (nrow(W) != nrow(V))
    stop("de novo profiles have ", nrow(W), " rows, expected ", nrow(V))
  if (!all(is.finite(W))) stop("de novo profiles contain non-finite values")

  # translate signeR's channel expression back to COSMIC, and verify
  if (!identical(rownames(W), rownames(V))) {
    conv <- signer_to_cosmic(rownames(W))
    if (!is.null(conv) && identical(conv, rownames(V))) {
      rownames(W) <- conv
    } else if (setequal(rownames(W), rownames(V))) {
      W <- W[rownames(V), , drop = FALSE]
    } else {
      stop("de novo profiles are not on the same 96 channels as V")
    }
  }
  if (is.null(colnames(W))) colnames(W) <- paste0("D", seq_len(ncol(W)))
  W
}

# MappingSignature() returns the selected catalogue names in $Reference
mapped_signatures <- function(W_hat, W_ref) {
  m <- SATS::MappingSignature(W_hat = W_hat, W_ref = W_ref)
  if (is.null(m) || !NROW(m)) return(list(names = character(0), table = m))
  nm <- if ("Reference" %in% names(m)) as.character(m$Reference) else {
    hit <- character(0)
    for (cn in setdiff(names(m), "freq")) {
      v <- as.character(m[[cn]])
      if (any(v %in% colnames(W_ref))) { hit <- v; break }
    }
    hit
  }
  nm <- unique(intersect(nm, colnames(W_ref)))
  list(names = nm, table = m)
}

## ------------------------------------------------------------------- output --

# SATS ships SBS2 and SBS13 merged as SBS2_13 because their activities are
# collinear. Split the burden evenly so the output carries the same signature
# names as every other tool.
split_apobec <- function(burden) {
  if (!"SBS2_13" %in% rownames(burden)) return(burden)
  half <- burden["SBS2_13", , drop = FALSE] / 2
  out <- burden[rownames(burden) != "SBS2_13", , drop = FALSE]
  out <- rbind(out, SBS2 = half[1, ], SBS13 = half[1, ])
  out[order(rownames(out)), , drop = FALSE]
}

## ------------------------------------------------------------------- driver --

run_sats_by_cohort <- function(
    maf, results_dir, genomic_info, clinical,
    cohort_column       = "CANCER_TYPE",
    ranks               = 1:12,
    pool_size           = 100L,
    min_pooled_profiles = 30L,
    min_cohort_samples  = 30L,
    min_signatures      = 2L,
    n_start             = 50L,
    seed                = 20260828L,
    testing_burn        = 500L,
    testing_eval        = 500L,
    main_burn           = 3000L,
    main_eval           = 1000L,
    parallelization     = "multisession",
    resume              = TRUE,
    out_name            = "SATS_ByCohort_Cosmic_V3.4") {

  out_dir <- file.path(results_dir, out_name)
  ckpt    <- file.path(out_dir, "checkpoints")
  dir.create(ckpt, recursive = TRUE, showWarnings = FALSE)
  channels <- sbs96_order()

  ## V and L, once for every tumour ------------------------------------------
  message("reading ", maf)
  mt <- data.table::fread(maf, sep = "\t", quote = "", header = TRUE,
                          showProgress = FALSE)
  need <- c("Chromosome", "Start_Position", "Reference_Allele",
            "Tumor_Seq_Allele2", "Variant_Type", "Tumor_Sample_Barcode")
  if (!all(need %in% names(mt)))
    stop("MAF must contain: ", paste(setdiff(need, names(mt)), collapse = ", "))
  message("  ", nrow(mt), " rows, ", uniqueN(mt$Tumor_Sample_Barcode), " tumours")

  V <- build_V(mt, channels, samples = unique(as.character(mt$Tumor_Sample_Barcode)))
  L <- build_L(genomic_info, clinical, colnames(V), channels,
               cache = file.path(out_dir, "panel_context.rds"))

  # keep only tumours present in both
  both <- intersect(colnames(V), colnames(L))
  if (!length(both)) stop("no tumour has both counts and an opportunity column")
  if (length(both) < ncol(V))
    message("  ", ncol(V) - length(both), " tumours dropped: no assay in the ",
            "clinical table or no opportunity for their assay")
  V <- V[, both, drop = FALSE]; L <- L[, both, drop = FALSE]
  stopifnot(identical(rownames(V), channels), identical(rownames(L), channels),
            identical(colnames(V), colnames(L)))
  message("  V and L: ", nrow(V), " x ", ncol(V))

  ## cohorts -----------------------------------------------------------------
  cl <- data.table::fread(clinical, sep = "\t", quote = "", skip = "PATIENT_ID",
                          showProgress = FALSE)
  if (!cohort_column %in% names(cl))
    stop("clinical table has no column '", cohort_column, "'")
  ct <- setNames(as.character(cl[[cohort_column]]), cl$SAMPLE_ID)[colnames(V)]
  ct[is.na(ct) | ct == ""] <- NA_character_
  if (anyNA(ct))
    message("  ", sum(is.na(ct)), " tumours have no ", cohort_column,
            " and are excluded")
  groups <- split(colnames(V)[!is.na(ct)], ct[!is.na(ct)])
  groups <- groups[order(-lengths(groups))]
  message("  cohorts: ", length(groups))

  e <- new.env(); utils::data("RefTMB", package = "SATS", envir = e)
  Wref <- as.matrix(get("RefTMB", envir = e)$TMB_SBS_v3.4)[channels, , drop = FALSE]

  ## per cohort --------------------------------------------------------------
  summ <- list(); siglong <- list(); burdens <- list()
  for (cname in names(groups)) {
    cols <- groups[[cname]]
    slug <- gsub("[^A-Za-z0-9]+", "_", cname)
    cfile <- file.path(ckpt, paste0(slug, ".rds"))

    if (resume && file.exists(cfile)) {
      got <- readRDS(cfile)
      message("\n--- ", cname, ": restored from checkpoint (", got$rec$status, ")")
      summ[[cname]] <- got$rec
      if (!is.null(got$B)) burdens[[cname]] <- got$B
      if (!is.null(got$sigs) && length(got$sigs))
        siglong[[cname]] <- data.table(cohort = cname, signature = got$sigs)
      next
    }

    n <- length(cols)
    rec <- data.table(cohort = cname, n_tumours = n, n_fitted = NA_integer_,
                      n_pools = NA_integer_, K_de_novo = NA_integer_,
                      n_signatures = NA_integer_, signatures = NA_character_,
                      status = NA_character_)
    message("\n--- ", cname, ": ", n, " tumours")
    if (n < min_cohort_samples) {
      rec$status <- sprintf("skipped (< min_cohort_samples %d)", min_cohort_samples)
      message("  ", rec$status)
      summ[[cname]] <- rec; saveRDS(list(rec = rec), cfile); next
    }

    Vc <- V[, cols, drop = FALSE]; Lc <- L[, cols, drop = FALSE]
    keep <- colSums(Vc) >= 1L
    rec$n_fitted <- sum(keep)
    if (!any(keep)) {
      rec$status <- "no tumour with a usable SNV"
      summ[[cname]] <- rec; saveRDS(list(rec = rec), cfile); next
    }
    Vf <- Vc[, keep, drop = FALSE]; Lf <- Lc[, keep, drop = FALSE]

    # pool, shrinking the group size rather than dropping below the floor
    m <- ncol(Vf); ps <- pool_size
    if (m %/% ps < min_pooled_profiles) ps <- max(1L, m %/% min_pooled_profiles)
    if (ps <= 1L) {
      Vd <- Vf; Ld <- Lf
      message("  discovery on ", ncol(Vd), " individual tumours")
    } else {
      pid <- ceiling(seq_len(m) / ps)
      Vd <- t(rowsum(t(Vf), group = pid, reorder = FALSE))
      Ld <- t(rowsum(t(Lf), group = pid, reorder = FALSE))
      message("  discovery on ", ncol(Vd), " pooled profiles (groups of ", ps, ")")
    }
    rec$n_pools <- ncol(Vd)

    W_hat <- tryCatch(
      de_novo(Vd, Ld, ranks = ranks, seed = seed,
              testing_burn = testing_burn, testing_eval = testing_eval,
              main_burn = main_burn, main_eval = main_eval,
              parallelization = parallelization),
      error = function(err) { message("  de novo FAILED: ", conditionMessage(err)); NULL })
    if (is.null(W_hat)) {
      rec$status <- "de novo failed"
      summ[[cname]] <- rec; saveRDS(list(rec = rec), cfile); next
    }
    rec$K_de_novo <- ncol(W_hat)

    mp <- tryCatch(mapped_signatures(W_hat, Wref),
                   error = function(err) {
                     message("  mapping FAILED: ", conditionMessage(err))
                     list(names = character(0), table = NULL) })
    SBS.list <- mp$names
    rec$n_signatures <- length(SBS.list)
    rec$signatures <- paste(SBS.list, collapse = ",")
    if (length(SBS.list) < min_signatures) {
      rec$status <- sprintf("degenerate (%d < min_signatures %d)",
                            length(SBS.list), min_signatures)
      message("  ", rec$status)
      summ[[cname]] <- rec; saveRDS(list(rec = rec), cfile); next
    }
    message("  K = ", ncol(W_hat), " -> ", length(SBS.list), ": ",
            paste(SBS.list, collapse = ", "))

    W_star <- Wref[, SBS.list, drop = FALSE]
    fit <- tryCatch(
      SATS::EstimateSigActivity(V = Vf, L = Lf, W = W_star, n.start = n_start),
      error = function(err) { message("  refit FAILED: ", conditionMessage(err)); NULL })
    if (is.null(fit)) {
      rec$status <- "refit failed"
      summ[[cname]] <- rec; saveRDS(list(rec = rec), cfile); next
    }
    bf <- SATS::CalculateSignatureBurdens(L = Lf, W = W_star, H = fit$H)

    B <- matrix(0, ncol(Wref), length(cols), dimnames = list(colnames(Wref), cols))
    B[rownames(bf), colnames(bf)] <- bf
    rec$status <- "ok"
    burdens[[cname]] <- B
    summ[[cname]] <- rec
    siglong[[cname]] <- data.table(cohort = cname, signature = SBS.list)
    saveRDS(list(rec = rec, B = B, sigs = SBS.list, mapped = mp$table), cfile)
    if (!is.null(mp$table))
      data.table::fwrite(as.data.frame(mp$table),
                         file.path(out_dir, paste0("mapped_", slug, ".tsv")), sep = "\t")
  }

  ## combine -----------------------------------------------------------------
  S <- data.table::rbindlist(summ, fill = TRUE)
  data.table::fwrite(S, file.path(out_dir, "cohort_summary.tsv"), sep = "\t")
  if (!length(burdens)) {
    cat("\nno cohort produced a fit\n"); print(S, row.names = FALSE)
    return(invisible(list(summary = S)))
  }
  cols_all <- unlist(lapply(burdens, colnames), use.names = FALSE)
  burden <- matrix(0, ncol(Wref), length(cols_all),
                   dimnames = list(colnames(Wref), cols_all))
  for (B in burdens) burden[, colnames(B)] <- B
  burden <- split_apobec(burden)

  out <- data.table::data.table(sample_name = colnames(burden))
  out <- cbind(out, data.table::as.data.table(t(burden)))
  data.table::fwrite(out, file.path(out_dir, "SATS_assignment_by_cohort.csv"),
                     sep = ",", row.names = FALSE)
  SL <- data.table::rbindlist(siglong, fill = TRUE)
  if (nrow(SL))
    data.table::fwrite(SL, file.path(out_dir, "cohort_signatures.tsv"), sep = "\t")

  cat("\n================ COHORT SUMMARY ================\n")
  print(S, row.names = FALSE)
  cat("\ncohorts fitted: ", sum(S$status == "ok", na.rm = TRUE), " of ", nrow(S),
      "\ntumours with burdens: ", ncol(burden),
      "\nwritten: ", file.path(out_dir, "SATS_assignment_by_cohort.csv"), "\n", sep = "")
  invisible(list(out = out, summary = S, signatures = SL))
}


## ------------------------------------------------------- standard entry point

# The contract the other assignment scripts follow:
#
#   run_<tool>(maf, results_dir, ...)  writes
#   <results_dir>/<Tool>_Cosmic_V3.4/<Tool>_assignment.csv
#
# one row per tumour, one column per COSMIC v3.4 signature, values are burdens.
# SATS additionally needs the target territory and the sample-to-assay map --
# without them there is no opportunity matrix and no SATS -- so genomic_info and
# clinical are required rather than optional.
run_sats <- function(maf, results_dir, genomic_info = NULL, clinical = NULL, ...) {
  if (is.null(genomic_info) || is.null(clinical))
    stop("run_sats() needs `genomic_info` (target intervals) and `clinical` ",
         "(SAMPLE_ID + SEQ_ASSAY_ID). SATS models the panel opportunity ",
         "explicitly; there is no sensible default for it.")

  res <- run_sats_by_cohort(maf = maf, results_dir = results_dir,
                            genomic_info = genomic_info, clinical = clinical, ...)
  if (is.null(res$out)) {
    warning("no cohort produced a fit; no assignment file written")
    return(invisible(res))
  }

  out_dir <- file.path(results_dir, "SATS_Cosmic_V3.4")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(out_dir, "SATS_assignment.csv")
  data.table::fwrite(res$out, path, sep = ",", row.names = FALSE)
  message("written: ", path)
  invisible(res)
}


run_sats(
  maf          = "../../../data/input/19.0_public/data_mutations_extended.txt",
  results_dir  = "test_sats/",
  genomic_info = "../../../data/input/19.0_public/genomic_information.txt",
  clinical     = "../../../data/input/19.0_public/data_clinical_sample.txt",
  cohort_column = "CANCER_TYPE")
