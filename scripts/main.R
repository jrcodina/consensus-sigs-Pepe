#!/usr/bin/env Rscript

library(optparse)
source("deconstructsigs_assignment.R")
source("sigminer_sigprofiler_assignment.R")
source("mesica_assignment.R")
source("sats_assignment.R")
source("consensus_pipeline.R")

option_list <- list(
  make_option(c("-m", "--maf"), type = "character", help = "Path to MAF input file", metavar = "FILE"),
  make_option(c("-o", "--output"), type = "character", help = "Output directory", metavar = "DIR"),
  make_option(c("-g", "--genomic-info"), type = "character", default = NULL, dest = "genomic_info",
              help = "Target intervals for the assay: a GENIE genomic_information.txt (1-based inclusive) or a plain 3-column BED (0-based). Used by SATS to build the mutation opportunity matrix.", metavar = "FILE"),
  make_option(c("-c", "--clinical"), type = "character", default = NULL,
              help = "Sample table with SAMPLE_ID and SEQ_ASSAY_ID, required only when --genomic-info covers more than one assay, so each sample gets its own panel's opportunity.", metavar = "FILE")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$maf) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("All options must be provided: --maf, --output")
}

results_dir <- file.path(opt$output, "consensus_sigs_results")
for (subdir in c("SigProfiler_Cosmic_V3.4", "SigMiner_Cosmic_V3.4", "deconstructSigs_Cosmic_V3.4", "MESiCA_Cosmic_V3.4", "SATS_Cosmic_V3.4", "processed", "consensus")) {
  dir.create(file.path(results_dir, subdir), recursive = TRUE, showWarnings = FALSE)
}

run_deconstructsigs(maf = opt$maf, results_dir = results_dir)
run_sigminer_sigprofiler(maf = opt$maf, results_dir = results_dir)
run_mesica(maf = opt$maf, results_dir = results_dir)
run_sats(maf = opt$maf, results_dir = results_dir, genomic_info = opt$genomic_info, clinical = opt$clinical)
run_consensus_pipeline(results_dir = results_dir)

cat("All signature assignments completed successfully.\n")
