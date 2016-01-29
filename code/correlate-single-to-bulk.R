#!/usr/bin/env Rscript

# Run -h for command-line options.

suppressMessages(library("docopt"))

"Correlate the expression in single cells to the bulk sample.

Usage:
correlate-single-cell-to-bulk.R [options] [--quantiles=<q>...] <num_cells> <seed> <single> <bulk>

Options:
  -h --help              Show this screen.
  --individual=<ind>     Only use data from ind, e.g. 19098
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

  # Load single cell data
  single_cells <- read.table(single_fname, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE)
  # Filter by individual
  if (!is.null(individual)) {
    single_cells <- single_cells[single_cells$individual == individual, ]
  }
  assert("Single cell data does not contain bulk samples",
         single_cells$well != "bulk")
  # Add rownames
  rownames(single_cells) <- paste(single_cells$individual, single_cells$replicate,
                                  single_cells$well, sep = ".")
  # Remove meta-info cols
  single_cells <- single_cells[, grepl("ENSG", colnames(single_cells)) |
                                 grepl("ERCC", colnames(single_cells))]
  # Transpose
  single_cells <- t(single_cells)
  # Fix ERCC names
  rownames(single_cells) <- sub(pattern = "\\.", replacement = "-",
                                rownames(single_cells))
  # Filter genes
  if (!is.null(keep_genes)) {
    keep_genes_list <- scan(keep_genes, what = "character", quiet = TRUE)
    single_cells <- single_cells[rownames(single_cells) %in% keep_genes_list, ]
  }
  # Keep only good quality cells
  if (!is.null(good_cells)) {
    assert("File with list of good quality cells exists.",
           file.exists(good_cells))
    good_cells_list <- scan(good_cells, what = "character", quiet = TRUE)
    single_cells <- single_cells[, colnames(single_cells) %in% good_cells_list]
    assert("There are quality cells to perform the analysis.",
           ncol(single_cells) > 0)
  }

  # Subsample number of single cells
  if (ncol(single_cells) < num_cells) {
    cat(sprintf("%d\t%d\tNA\tNA\tNA\n", num_cells, seed))
    return(invisible())
  }
  set.seed(seed)
  single_cells <- single_cells[, sample(1:ncol(single_cells), size = num_cells)]
  # Calculate cpm
  single_cells <- cpm(single_cells)
#   single_cells[1:10, 1:10]
#   dim(single_cells)

  # Load bulk cell data
  bulk_cells <- read.table(bulk_fname, header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE)
  # Filter by individual
  if (!is.null(individual)) {
    bulk_cells <- bulk_cells[bulk_cells$individual == individual, ]
  }
  assert("Bulk data does not contain single cell samples",
         bulk_cells$well == "bulk")
  # Add rownames
  rownames(bulk_cells) <- paste(bulk_cells$individual, bulk_cells$replicate,
                                  bulk_cells$well, sep = ".")
  # Remove meta-info cols
  bulk_cells <- bulk_cells[, grepl("ENSG", colnames(bulk_cells)) |
                             grepl("ERCC", colnames(bulk_cells))]
  # Transpose
  bulk_cells <- t(bulk_cells)
  # Fix ERCC names
  rownames(bulk_cells) <- sub(pattern = "\\.", replacement = "-",
                                rownames(bulk_cells))
  # Filter genes
  if (!is.null(keep_genes)) {
    bulk_cells <- bulk_cells[rownames(bulk_cells) %in% keep_genes_list, ]
  }
  # Calculate cpm
  bulk_cells <- cpm(bulk_cells)

#   bulk_cells[1:10, 1:10]
#   dim(bulk_cells)
#   head(bulk_cells)

  assert("Same number of genes in bulk and single cells.",
         nrow(bulk_cells) == nrow(single_cells))
  assert("Same order of genes in bulk and single cells.",
         rownames(bulk_cells) == rownames(single_cells))

  mean_single_expression <- log2(rowMeans(single_cells) + 1)
  mean_bulk_expression <- log2(rowMeans(bulk_cells) + 1)

  quantiles <- c(quantiles, 1)
  quantiles <- sort(quantiles)
  q_cutoffs <- quantile(mean_bulk_expression, probs = quantiles)
  q_r <- numeric(length = length(quantiles))
  q_n <- numeric(length = length(quantiles))
  for (i in 1:length(quantiles)) {
    # Correlate
    gene_in_quantile <- mean_bulk_expression <= q_cutoffs[i]
    q_n[i] <- sum(gene_in_quantile)
    q_r[i] <- cor(mean_single_expression[gene_in_quantile],
                  mean_bulk_expression[gene_in_quantile])
    # Output
    cat(sprintf("%d\t%d\t%f\t%f\t%d\n", num_cells, seed, quantiles[i], q_r[i], q_n[i]))
    # Remove genes already analyzed
    mean_single_expression <- mean_single_expression[!gene_in_quantile]
    mean_bulk_expression <- mean_bulk_expression[!gene_in_quantile]
  }
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
