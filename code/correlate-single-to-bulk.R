#!/usr/bin/env Rscript

# Run -h for command-line options.

suppressMessages(library("docopt"))

"Correlate the expression in single cells to the bulk sample.

Usage:
correlate-single-cell-to-bulk.R [options] [--quantiles=<q>...] <num_cells> <seed> <single> <bulk>

Options:
  -h --help              Show this screen.
  --individual=<ind>     Only use data from ind, e.g. NA19098
  --good_cells=<file>    A 1-column file with the names of good quality cells to maintain
  --keep_genes=<file>    A 1-column file with the names of genes to maintain
  -q --quantiles=<q>     Calculate the correlation for the genes separated by the provided
                         quantiles, e.g. -q .25 -q .75

Arguments:
  num_cells     number of single cells to subsample
  seed          seed for random number generator
  single        sample-by-gene matrix of single cell data
  bulk          sample-by-gene matrix of bulk cell data" -> doc

main <- function(num_cells, seed, single_fname, bulk_fname, individual = NULL,
                 good_cells = NULL, keep_genes = NULL, quantiles = NULL) {
  suppressPackageStartupMessages(library("edgeR"))
  library("testit")
  id <- "single-to-bulk-correlation"

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

  # Load single cell data
  single_cells <- read.table(single_fname, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE)
  assert("Single cell data does not contain bulk samples",
         single_cells$well != "bulk")

  # Filter individuals, cells, and genes. Also transpose to gene-by-sample.
  single_cells <- prepare_counts(single_cells, individual = individual,
                                 good_cells_list = good_cells_list,
                                 keep_genes_list = keep_genes_list)
  # Subsample number of single cells
  if (ncol(single_cells) < num_cells) {
    cat(sprintf("%d\t%d\tNA\tNA\tNA\n", num_cells, seed))
    return(invisible())
  }
  set.seed(seed)
  single_cells <- single_cells[, sample(1:ncol(single_cells), size = num_cells)]

#   single_cells[1:10, 1:10]
#   dim(single_cells)

  # Load bulk cell data
  bulk_cells <- read.table(bulk_fname, header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE)
  assert("Bulk data does not contain single cell samples",
         bulk_cells$well == "bulk")
  # Filter individuals and genes. Also transpose to gene-by-sample.
  bulk_cells <- prepare_counts(bulk_cells, individual = individual,
                               keep_genes_list = keep_genes_list)

#   bulk_cells[1:10, 1:10]
#   dim(bulk_cells)
#   head(bulk_cells)

  assert("Same number of genes in bulk and single cells.",
         nrow(bulk_cells) == nrow(single_cells))
  assert("Same order of genes in bulk and single cells.",
         rownames(bulk_cells) == rownames(single_cells))

  # Do not include ERCC control genes
  endogenous <- grep("ENSG", rownames(single_cells))
  single_cells <- single_cells[endogenous, ]
  bulk_cells <- bulk_cells[endogenous, ]

  # For single cells, sum the counts across the single cells and then calculate
  # log2 cpm
  single_cells_sum <- as.data.frame(rowSums(single_cells))
  single_cells_sum_cpm <- cpm(single_cells_sum, log = TRUE, prior.count = 1)
  single_cells_sum_cpm <- as.numeric(single_cells_sum_cpm)
  # For bulk samples, calculate cpm for each replicate and then calculate the
  # mean across the replicates
  bulk_cells_cpm <- cpm(bulk_cells, log = TRUE, prior.count = 1)
  bulk_cells_cpm_mean <- rowMeans(bulk_cells_cpm)

  quantiles <- c(quantiles, 1)
  quantiles <- sort(quantiles)
  q_cutoffs <- quantile(bulk_cells_cpm_mean, probs = quantiles)
  q_r <- numeric(length = length(quantiles))
  q_n <- numeric(length = length(quantiles))
  for (i in 1:length(quantiles)) {
    # Correlate
    gene_in_quantile <- bulk_cells_cpm_mean <= q_cutoffs[i]
    q_n[i] <- sum(gene_in_quantile)
    q_r[i] <- cor(single_cells_sum_cpm[gene_in_quantile],
                  bulk_cells_cpm_mean[gene_in_quantile])
    # Output
    cat(sprintf("%d\t%d\t%f\t%f\t%d\n", num_cells, seed, quantiles[i], q_r[i], q_n[i]))
    # Remove genes already analyzed
    single_cells_sum_cpm <- single_cells_sum_cpm[!gene_in_quantile]
    bulk_cells_cpm_mean <- bulk_cells_cpm_mean[!gene_in_quantile]
  }
}

# Converts a sample-by-gene data frame to a filtered gene-by-sample data frame.
#
# x - sample-by-gene data frame
# individual - character vector of individuals to keep, e.g. NA19098
# good_cells_list - A character vector with the names of good quality cells to maintain
# keep_genes_list - A character vector with the names of genes to maintain
#
prepare_counts <- function(x, individual = NULL, good_cells_list = NULL,
                           keep_genes_list = NULL) {
  library("testit")
  assert("Input is a data frame", class(x) == "data.frame")
  assert("Input is not empty", dim(x) > 0)
  assert("Input has necessary columns",
         c("individual", "replicate", "well") %in% colnames(x))
  # Filter by individual
  if (!is.null(individual)) {
    x <- x[x$individual == individual, ]
  }
  # Add rownames
  rownames(x) <- paste(x$individual, x$replicate, x$well, sep = ".")
  # Remove meta-info cols
  x <- x[, grepl("ENSG", colnames(x)) | grepl("ERCC", colnames(x))]
  # Transpose
  x <- t(x)
  # Fix ERCC names
  rownames(x) <- sub(pattern = "\\.", replacement = "-", rownames(x))
  # Filter genes
  if (!is.null(keep_genes_list)) {
    x <- x[rownames(x) %in% keep_genes_list, ]
  }
  # Keep only good quality cells
  if (!is.null(good_cells_list)) {
    x <- x[, colnames(x) %in% good_cells_list]
    assert("There are quality cells to perform the analysis.",
           ncol(x) > 0)
  }
  return(x)
}

if (!interactive() & getOption('run.main', default = TRUE)) {
  opts <- docopt(doc)
  main(num_cells = as.numeric(opts$num_cells),
       seed = as.numeric(opts$seed),
       single_fname = opts$single,
       bulk_fname = opts$bulk,
       individual = opts$individual,
       good_cells = opts$good_cells,
       keep_genes = opts$keep_genes,
       quantiles = as.numeric(opts$quantiles))
} else if (interactive() & getOption('run.main', default = TRUE)) {
  # what to do if interactively testing
  main(num_cells = 20,
       seed = 1,
       single_fname = "/mnt/gluster/home/jdblischak/ssd/subsampled/counts-matrix/250000-molecules-raw-single-per-sample.txt",
       bulk_fname = "../data/reads-raw-bulk-per-sample.txt",
       individual = "NA19098",
       good_cells = "../data/quality-single-cells.txt",
       keep_genes = "../data/genes-pass-filter.txt",
       quantiles = c(.25, .5, .75))
}
