#!/usr/bin/env Rscript

# Run -h for command-line options.

library("testit")
suppressPackageStartupMessages(library("docopt"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("edgeR"))

"Calculate statistics on the subsampled data.

Usage:
subsampler.R [options] [--quantiles=<q>...] <num_cells> <seed> <single_sub> <single_full> <ercc> <bulk>

Options:
  -h --help              Show this screen.
  --individual=<ind>     Only use data from ind, e.g. NA19098
  --replicate=<rep>      Only use data from rep, e.g. r1
  --good_cells=<file>    A 1-column file with the names of good quality cells to maintain
  --keep_genes=<file>    A 1-column file with the names of genes to maintain
  -q --quantiles=<q>     Calculate the correlation for the genes separated by the provided
                         quantiles, e.g. -q .25 -q .75
  -o --outfile=<file>    Output file to write results (default: stdout)
  -d --diagnose          Saves diagnostic plots to PDF file. Name of file is derived
                         from outfile option, otherwise written to Rplots.pdf

Arguments:
  num_cells     number of single cells to subsample
  seed          seed for random number generator
  single_sub    sample-by-gene matrix of subsampled single cell data
  single_full   sample-by-gene matrix of full single cell data
  ercc          the expected ERCC counts for each single cell well
                  (column id contains ERCC identifiers and column
                  ercc_expected_well contains number of expected molecules)
  bulk          sample-by-gene matrix of bulk cell data" -> doc

# Main function ----------------------------------------------------------------
main <- function(num_cells, seed, single_sub_fname, single_full_fname,
                 ercc_fname, bulk_fname,
                 individual = NULL, replicate = NULL,
                 good_cells = NULL, keep_genes = NULL, quantiles = NULL,
                 outfile = NULL, diagnose = FALSE) {

  # Data frame to contain all results
  results <- data.frame(seed, subsampled_cells = num_cells)
  if (!is.null(individual)) {
    results$individual = individual
  } else {
    results$individual = NA
  }
  if (!is.null(replicate)) {
    results$replicate = replicate
  } else {
    results$replicate = NA
  }

  # If an output filename is specified, save the PDF using the same name with .pdf
  # added as the file extension.
  if (diagnose & !is.null(outfile)) {
    pdf_file <- paste0(outfile, ".pdf")
    pdf(pdf_file)
    on.exit(dev.off())
  }
  # If no output file is specified, send results to standard out
  if (is.null(outfile)) {
    outfile <- ""
  }

  # Load filtering data: good_cells and keep_genes
  if (!is.null(keep_genes)) {
    assert("File with list of genes to keep exists.",
           file.exists(keep_genes))
    keep_genes_list <- scan(keep_genes, what = "character", quiet = TRUE)
  } else {
    keep_genes_list <- NULL
  }
  if (!is.null(good_cells)) {
    assert("File with list of good quality cells exists.",
           file.exists(good_cells))
    good_cells_list <- scan(good_cells, what = "character", quiet = TRUE)
  } else {
    good_cells_list <- NULL
  }

  # Load count data sets. Filter individuals, cells, and genes. Also
  # transpose to gene-by-sample.

  # Subsampled single cell data
  single_cells_sub <- prepare_counts(single_sub_fname, individual = individual,
                                 replicate = replicate,
                                 good_cells_list = good_cells_list,
                                 keep_genes_list = keep_genes_list)
  # Full single cell data
  single_cells_full <- prepare_counts(single_full_fname, individual = individual,
                                      replicate = replicate,
                                      good_cells_list = good_cells_list,
                                      keep_genes_list = keep_genes_list)
  # bulk cell data
  bulk_cells <- prepare_counts(bulk_fname, individual = individual,
                               replicate = replicate,
                               keep_genes_list = keep_genes_list)

  assert("Same number of genes in bulk and single cells.",
         nrow(bulk_cells) == nrow(single_cells_sub),
         nrow(single_cells_full) == nrow(single_cells_sub))
  assert("Same order of genes in bulk and single cells.",
         rownames(bulk_cells) == rownames(single_cells_sub),
         rownames(single_cells_full) == rownames(single_cells_sub))
  assert("Subsampled cells are a subset of all single cells",
         colnames(single_cells_sub) %in% colnames(single_cells_full))

  # Some diagnostic plots exploring the total counts for endogenous and ERCC
  # genes between the subsampled and original full single cell data sets.
  if (diagnose) {
    totals_endo_single_sub <- colSums(single_cells_sub[
                              grep("ENSG", rownames(single_cells_sub)), ])
    totals_ercc_single_sub <- colSums(single_cells_sub[
                              grep("ERCC", rownames(single_cells_sub)), ])
    totals_endo_single_full <- colSums(single_cells_full[
                               grep("ENSG", rownames(single_cells_full)), ])
    totals_ercc_single_full <- colSums(single_cells_full[
                               grep("ERCC", rownames(single_cells_full)), ])
    op <- par(mfrow = c(2, 2))
    hist(totals_endo_single_sub, main = "Subsampled Endogenous")
    hist(totals_ercc_single_sub, main = "Subsampled ERCC")
    hist(totals_endo_single_full, main = "Full Endogenous")
    hist(totals_ercc_single_full, main = "Full ERCC")
    par(op)
    # Compare totals from full and subsampled data sets. Have to subset to only
    # include cells in the subsampled set.
    cells_in_sub <- colnames(single_cells_full) %in% colnames(single_cells_sub)
    totals_combined <- data.frame(endo_sub = totals_endo_single_sub,
                                  ercc_sub = totals_ercc_single_sub,
                                  endo_full = totals_endo_single_full[cells_in_sub],
                                  ercc_full = totals_ercc_single_full[cells_in_sub])
    pairs(totals_combined)
  }

  # potential cells - the number of single cells available after filtering for
  # individual, replicate, and quality
  results$potential_cells <- ncol(single_cells_full)
  # available cells - the number of potential cells that had sufficient
  # sequencing depth to be included in this subsample
  results$available_cells <- ncol(single_cells_sub)

  # Subsample number of single cells
  if (results$available_cells < num_cells) {
    write.table(results, file = outfile, quote = FALSE, row.names = FALSE,
                sep = "\t")
    return(invisible())
  }
  set.seed(seed)
  single_cells_sub <- single_cells_sub[, sample(1:ncol(single_cells_sub), size = num_cells)]

  # Split endogenous and ERCC genes
  ensg_index <- grepl("ENSG", rownames(bulk_cells))
  ercc_index <- grepl("ERCC", rownames(bulk_cells))
  assert("Both endogenous Ensembl and ERCC genes available for analysis",
         sum(ensg_index) > 0, sum(ercc_index) > 0)

  # Calculate mean gene expression:
  # For single cells, sum the counts across the single cells and then calculate
  # log2 cpm
  single_cells_cpm_mean_ensg <- calc_mean_single_cell(single_cells_sub[ensg_index, , drop = FALSE])
  single_cells_cpm_mean_ercc <- calc_mean_single_cell(single_cells_sub[ercc_index, , drop = FALSE])
  # For bulk samples, calculate cpm for each replicate and then calculate the
  # mean across the replicates
  bulk_cells_cpm_mean_ensg <- calc_mean_bulk_cell(bulk_cells[ensg_index, , drop = FALSE])
  bulk_cells_cpm_mean_ercc <- calc_mean_bulk_cell(bulk_cells[ercc_index, , drop = FALSE])

  # Calculate correlation between single cells and bulk samples
  results$mean_cor_ensg <- calc_mean_cor(single_cells_cpm_mean_ensg, bulk_cells_cpm_mean_ensg,
                                         diagnose = diagnose, method = "spearman")
  results$mean_cor_ercc <- calc_mean_cor(single_cells_cpm_mean_ercc, bulk_cells_cpm_mean_ercc,
                                         diagnose = diagnose, method = "spearman")

  write.table(results, file = outfile, quote = FALSE, row.names = FALSE,
              sep = "\t")
}

# Data import functions --------------------------------------------------------

# Converts a sample-by-gene data frame to a filtered gene-by-sample matrix.
#
# fname - filename of sample-by-gene count matrix
# individual - character vector of individuals to keep, e.g. NA19098
# good_cells_list - A character vector with the names of good quality cells to maintain
# keep_genes_list - A character vector with the names of genes to maintain
#
prepare_counts <- function(fname, individual = NULL, replicate = NULL,
                           good_cells_list = NULL, keep_genes_list = NULL) {
  assert("Input file exists", file.exists(fname))
  x <- fread(fname)
  setDF(x)
  assert("Input has necessary columns",
         c("individual", "replicate", "well") %in% colnames(x))
  # Filter by individual
  if (!is.null(individual)) {
    x <- x[x$individual == individual, ]
  }
  # Filter by replicate
  if (!is.null(replicate)) {
    x <- x[x$replicate == replicate, ]
  }
  # Add rownames
  rownames(x) <- paste(x$individual, x$replicate, x$well, sep = ".")
  # Remove meta-info cols
  x <- x[, grepl("ENSG", colnames(x)) | grepl("ERCC", colnames(x)), drop = FALSE]
  # Transpose
  x <- t(x)
  # Filter genes
  if (!is.null(keep_genes_list)) {
    x <- x[rownames(x) %in% keep_genes_list, , drop = FALSE]
  }
  # Keep only good quality cells
  if (!is.null(good_cells_list)) {
    x <- x[, colnames(x) %in% good_cells_list]
    assert("There are quality cells to perform the analysis.",
           ncol(x) > 0)
  }
  assert("Output is a matrix", class(x) == "matrix")
  return(x)
}

# Utility functions ------------------------------------------------------------

# Calculate mean expression per gene across single cells.
#
# Method: First sum counts across single cells and then calcalute log2 cpm using
# edgeR::cpm.
#
# x - a gene-by-sample matrix of gene counts
#
# Returns a numeric vector of mean counts per million
#
calc_mean_single_cell <- function(x) {
  assert("Input is matrix", is.matrix(x))
  assert("Input is not empty", dim(x) > 0)
  single_cells_sum <- as.data.frame(rowSums(x))
  single_cells_sum_cpm <- cpm(single_cells_sum, log = TRUE)
  single_cells_sum_cpm <- as.numeric(single_cells_sum_cpm)
  assert("No strange data values in output", !is.na(single_cells_sum_cpm),
         !is.null(single_cells_sum_cpm), !is.nan(single_cells_sum_cpm),
         !is.infinite(single_cells_sum_cpm))
  return(single_cells_sum_cpm)
}

# Calculate mean expression per gene across bulk cells.
#
# Method: First calculate log2 cpm (edgeR::cpm) for each replicate and then calculate the mean across
# the replicates.
#
# x - a gene-by-sample matrix of gene counts
#
# Returns a numeric vector of mean counts per million
#
calc_mean_bulk_cell <- function(x) {
  assert("Input is matrix", is.matrix(x))
  assert("Input is not empty", dim(x) > 0)
  bulk_cells_cpm <- cpm(x, log = TRUE)
  bulk_cells_cpm_mean <- rowMeans(bulk_cells_cpm)
  assert("No strange data values in output", !is.na(bulk_cells_cpm_mean),
         !is.null(bulk_cells_cpm_mean), !is.nan(bulk_cells_cpm_mean),
         !is.infinite(bulk_cells_cpm_mean))
  return(bulk_cells_cpm_mean)
}

# Calculate correlation between two numeric vectors.
#
# x - a numeric vector
# y - a numeric vector
# method - method to caclulate correlation, see ?cor for options.
# diagnose - Create diagnostic plot
# digits - The number of digits for rounding. Default is 5.
#
# Returns a numeric vector of mean counts per million
#
calc_mean_cor <- function(x, y, method = "pearson", diagnose = FALSE,
                          digits = 5) {
  assert("Proper input", length(x) == length(y), is.numeric(x), is.numeric(y))
  correlation <- cor(x, y, method = method)
  if (diagnose) {
    model <- lm(y ~ x)
    plot(x, y)
    abline(0, 1, col = "red")
    abline(model, col = "blue")
    title(main = sprintf("r: %.4f", correlation),
          sub = "Red: 1-1    Blue: Best fit")
  }
  assert("Proper output", correlation >= -1, correlation <= 1,
         length(correlation) == 1)
  correlation <- round(correlation, digits = digits)
  return(correlation)
}


# ------------------------------------------------------------------------------
# Input parameters

if (!interactive() & getOption('run.main', default = TRUE)) {
  opts <- docopt(doc)
  main(num_cells = as.numeric(opts$num_cells),
       seed = as.numeric(opts$seed),
       single_sub_fname = opts$single_sub,
       single_full_fname = opts$single_full,
       ercc = opts$ercc,
       bulk_fname = opts$bulk,
       individual = opts$individual,
       replicate = opts$replicate,
       good_cells = opts$good_cells,
       keep_genes = opts$keep_genes,
       quantiles = as.numeric(opts$quantiles),
       outfile = opts$outfile,
       diagnose = opts$diagnose)
} else if (interactive() & getOption('run.main', default = TRUE)) {
  # what to do if interactively testing
  main(num_cells = 20,
       seed = 1,
       single_sub_fname = "/mnt/gluster/home/jdblischak/ssd/subsampled/counts-matrix/250000-molecules-raw-single-per-sample.txt",
       single_full_fname = "../data/molecules-raw-single-per-sample.txt",
       ercc = "../data/expected-ercc-molecules.txt",
       bulk_fname = "../data/reads-raw-bulk-per-sample.txt",
       individual = "NA19098",
       replicate = "r1",
       good_cells = "../data/quality-single-cells.txt",
       keep_genes = "../data/genes-pass-filter.txt",
       quantiles = c(.25, .5, .75),
       diagnose = TRUE)
}
