# SATS assignment
#
# SATS (Zhu et al.) differs from the other four tools in one respect: it models
# the mutation *opportunity* explicitly.  Where deconstructSigs, sigminer and
# SigProfiler are handed a whole-exome correction, SATS takes an opportunity
# matrix L as a separate input and fits
#
#     V[p,n] ~ Poisson( L[p,n] * sum_k W[p,k] * H[k,n] )
#
# with W fixed (COSMIC v3.4, TMB-normalised, shipped with the package) and only
# the activities H estimated.  Burdens -- expected mutation counts per signature
# per sample -- are the output, the same kind of table the other assignment
# scripts produce.
#
# L is built from the assay's own target intervals, supplied with --genomic-info.
# Two formats are accepted:
#
#   * GENIE genomic_information.txt -- header with Chromosome, Start_Position,
#     End_Position and optionally SEQ_ASSAY_ID.  Coordinates are 1-based
#     inclusive, which is NOT BED convention, so no offset is applied.
#   * a plain BED -- three headerless columns, 0-based half-open, so the start
#     is shifted by one.
#
# If the interval file covers more than one assay, --clinical must also be given
# (a GENIE data_clinical_sample.txt, or any table with SAMPLE_ID and
# SEQ_ASSAY_ID) so each sample gets its own panel's opportunity.  Without an
# interval file the script falls back to the exome composition used by the other
# tools and warns, because on panel data that is a known bias, not a neutral
# default.
#
# Output: <results_dir>/SATS_Cosmic_V3.4/SATS_assignment.csv
#         one row per sample, one column per COSMIC v3.4 SBS signature.

library(data.table)
library(dplyr)
library(GenomicRanges)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(deconstructSigs)
library(SATS)

# ---------------------------------------------------------------- helpers ----

# the 96 SBS channels in COSMIC order, e.g. "A[C>A]A"
sats_sbs96_order <- function() {
  subs <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  flank <- c("A", "C", "G", "T")
  unlist(lapply(subs, function(s) {
    as.vector(t(outer(flank, flank, function(a, b) paste0(a, "[", s, "]", b))))
  }))
}

# the 32 pyrimidine trinucleotides, in the same order as the 96 channels
sats_sbs96_trinuc <- function(ch = sats_sbs96_order()) {
  paste0(substr(ch, 1, 1), substr(ch, 3, 3), substr(ch, 7, 7))
}

sats_revcomp <- function(x) as.character(reverseComplement(DNAStringSet(x)))

# 96 x N mutation matrix from a MAF, folded onto the pyrimidine strand
sats_build_V <- function(mut, genome, channels) {
  gr <- GRanges(
    seqnames = paste0("chr", mut$Chromosome),
    ranges = IRanges(start = mut$Start_Position - 1L, end = mut$Start_Position + 1L)
  )
  inb <- start(gr) >= 1L & end(gr) <= seqlengths(genome)[as.character(seqnames(gr))]
  gr <- gr[inb]
  mut <- mut[inb]

  ctx <- as.character(getSeq(genome, gr))

  # the middle base must agree with the reference allele the MAF reports
  agrees <- substr(ctx, 2, 2) == mut$Reference_Allele
  if (mean(agrees) < 0.99) {
    warning(sprintf(
      "only %.1f%% of SNVs match hg19 at the reference base -- check the genome build",
      100 * mean(agrees)
    ))
  }
  ctx <- ctx[agrees]
  mut <- mut[agrees]

  alt <- mut$Tumor_Seq_Allele2
  purine <- substr(ctx, 2, 2) %in% c("A", "G")
  ctx[purine] <- sats_revcomp(ctx[purine])
  alt[purine] <- chartr("ACGT", "TGCA", alt[purine])
  ref <- substr(ctx, 2, 2)

  mut[, channel := paste0(substr(ctx, 1, 1), "[", ref, ">", alt, "]", substr(ctx, 3, 3))]
  mut <- mut[channel %in% channels]

  tab <- mut[, .N, by = .(channel, Tumor_Sample_Barcode)]
  samples <- sort(unique(mut$Tumor_Sample_Barcode))
  V <- matrix(0L, length(channels), length(samples), dimnames = list(channels, samples))
  V[cbind(match(tab$channel, channels), match(tab$Tumor_Sample_Barcode, samples))] <- tab$N
  V
}

# read either a GENIE genomic_information.txt or a plain BED into one table
sats_read_intervals <- function(path) {
  iv <- data.table::fread(path, sep = "\t", header = TRUE, showProgress = FALSE)
  needed <- c("Chromosome", "Start_Position", "End_Position")

  if (!all(needed %in% names(iv))) {
    # plain BED: headerless, 0-based half-open, so the start shifts by one
    iv <- data.table::fread(
      path, sep = "\t", header = FALSE, select = 1:3,
      col.names = needed, showProgress = FALSE
    )
    iv[, Start_Position := Start_Position + 1L]
    message("--genomic-info read as a 3-column BED (0-based); start coordinates shifted by 1")
  } else {
    message("--genomic-info read as a GENIE interval table (1-based inclusive); no offset applied")
  }

  iv[, Chromosome := sub("^chr", "", as.character(Chromosome))]
  iv[Chromosome %in% c(as.character(1:22), "X", "Y")]
}

# opportunity per Mb for the 96 channels, from a set of target intervals
sats_interval_opportunity <- function(iv, channels, genome) {
  gr <- GRanges(
    seqnames = paste0("chr", iv$Chromosome),
    # widen by one base so contexts straddling a target edge are counted
    ranges = IRanges(start = pmax(iv$Start_Position - 1L, 1L), end = iv$End_Position + 1L)
  )
  # rows overlap in GENIE's table -- a base sequenced twice is still one opportunity
  gr <- GenomicRanges::trim(GenomicRanges::reduce(gr))

  freq <- colSums(trinucleotideFrequency(getSeq(genome, gr)))
  all_tri <- names(freq)
  pyr <- all_tri[substr(all_tri, 2, 2) %in% c("C", "T")]
  folded <- freq[pyr] + freq[sats_revcomp(pyr)]
  names(folded) <- pyr
  as.numeric(folded[sats_sbs96_trinuc(channels)]) / 1e6
}

# the exome fallback, matching what the other tools apply
sats_exome_opportunity <- function(channels) {
  counts <- deconstructSigs::tri.counts.exome
  opp <- setNames(as.numeric(counts[[1]]), rownames(counts))
  as.numeric(opp[sats_sbs96_trinuc(channels)]) / 1e6
}

# 96 x N opportunity matrix, one column per sample
sats_build_L <- function(V, channels, genome, genomic_info = NULL, clinical = NULL) {
  samples <- colnames(V)

  if (is.null(genomic_info)) {
    warning(
      "no --genomic-info supplied: falling back to exome opportunity. ",
      "On panel data this inflates CpG-driven signatures (SBS1 in particular), ",
      "because a gene panel is markedly CpG-richer than the exome."
    )
    col <- sats_exome_opportunity(channels)
    return(matrix(col, length(channels), length(samples),
                  dimnames = list(channels, samples)))
  }

  iv <- sats_read_intervals(genomic_info)
  assays <- if ("SEQ_ASSAY_ID" %in% names(iv)) sort(unique(iv$SEQ_ASSAY_ID)) else character(0)

  # one assay, or no assay column: a single opportunity vector for everyone
  if (length(assays) <= 1L) {
    col <- sats_interval_opportunity(iv, channels, genome)
    message(sprintf("opportunity from %s intervals covering %.3f Mb",
                    format(nrow(iv), big.mark = ","), sum(col) / 3))
    return(matrix(col, length(channels), length(samples),
                  dimnames = list(channels, samples)))
  }

  # several assays: each sample needs its own panel's opportunity
  if (is.null(clinical)) {
    stop(
      "--genomic-info covers ", length(assays), " assays (SEQ_ASSAY_ID), so each ",
      "sample needs its own opportunity. Supply --clinical with a table of ",
      "SAMPLE_ID and SEQ_ASSAY_ID, or filter the interval file to one assay. ",
      "Pooling intervals across assays would give every sample the wrong L."
    )
  }

  clin <- data.table::fread(clinical, sep = "\t", quote = "",
                            skip = "PATIENT_ID", showProgress = FALSE)
  if (!all(c("SAMPLE_ID", "SEQ_ASSAY_ID") %in% names(clin)))
    stop("--clinical must contain SAMPLE_ID and SEQ_ASSAY_ID columns")

  panel_of <- setNames(as.character(clin$SEQ_ASSAY_ID), clin$SAMPLE_ID)[samples]
  unmapped <- is.na(panel_of) | !(panel_of %in% assays)

  by_assay <- vapply(
    assays,
    function(a) sats_interval_opportunity(iv[SEQ_ASSAY_ID == a], channels, genome),
    numeric(length(channels))
  )

  if (any(unmapped)) {
    warning(sum(unmapped), " of ", length(samples), " samples have no assay in ",
            "--clinical; they get the pooled opportunity across all assays")
    pooled <- sats_interval_opportunity(iv, channels, genome)
    by_assay <- cbind(by_assay, `__pooled__` = pooled)
    panel_of[unmapped] <- "__pooled__"
  }

  message(sprintf("opportunity built per assay for %d assays across %d samples",
                  length(assays), length(samples)))
  L <- by_assay[, panel_of, drop = FALSE]
  dimnames(L) <- list(channels, samples)
  L
}

# SATS ships SBS2 and SBS13 merged as one column, because their activities are
# collinear; the merged profile is an even average of the two.  Split the burden
# evenly so the output carries the same signature names as every other tool.
sats_split_apobec <- function(burden) {
  if (!"SBS2_13" %in% rownames(burden)) return(burden)
  half <- burden["SBS2_13", , drop = FALSE] / 2
  out <- burden[rownames(burden) != "SBS2_13", , drop = FALSE]
  out <- rbind(out, SBS2 = half[1, ], SBS13 = half[1, ])
  out[order(rownames(out)), , drop = FALSE]
}

# -------------------------------------------------------------------- run ----

run_sats <- function(maf, results_dir, genomic_info = NULL, clinical = NULL,
                     n_start = 30L) {
  genome <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
  channels <- sats_sbs96_order()

  maf_table <- data.table::fread(maf, sep = "\t", quote = "", header = TRUE)

  mut <- maf_table[
    Variant_Type == "SNP" &
      Reference_Allele %in% c("A", "C", "G", "T") &
      Tumor_Seq_Allele2 %in% c("A", "C", "G", "T") &
      Reference_Allele != Tumor_Seq_Allele2 &
      as.character(Chromosome) %in% c(as.character(1:22), "X"),
    .(Chromosome = as.character(Chromosome), Start_Position,
      Reference_Allele, Tumor_Seq_Allele2, Tumor_Sample_Barcode)
  ]
  if (nrow(mut) == 0L) stop("no usable SNVs in ", maf)

  V <- sats_build_V(mut, genome, channels)
  L <- sats_build_L(V, channels, genome, genomic_info = genomic_info, clinical = clinical)

  utils::data("RefTMB", package = "SATS", envir = environment())
  W <- as.matrix(RefTMB$TMB_SBS_v3.4)[channels, , drop = FALSE]

  # a sample with no mutations has nothing to fit; hold it out and return zeros
  informative <- colSums(V) >= 1L
  valid <- SATS::ValidateSATSInputs(
    V = V[, informative, drop = FALSE],
    L = L[, informative, drop = FALSE],
    W = W
  )
  fit <- SATS::EstimateSigActivity(
    V = valid$V, L = valid$L, W = valid$W, n.start = n_start
  )
  burden_fit <- SATS::CalculateSignatureBurdens(L = valid$L, W = valid$W, H = fit$H)

  burden <- matrix(
    0, nrow(burden_fit), ncol(V),
    dimnames = list(rownames(burden_fit), colnames(V))
  )
  burden[, colnames(burden_fit)] <- burden_fit
  burden <- sats_split_apobec(burden)

  out <- data.table::data.table(sample_name = colnames(burden))
  out <- cbind(out, data.table::as.data.table(t(burden)))

  dir.create(file.path(results_dir, "SATS_Cosmic_V3.4"), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(
    x = out, sep = ",", row.names = FALSE,
    file = file.path(results_dir, "SATS_Cosmic_V3.4", "SATS_assignment.csv")
  )

  invisible(out)
}
