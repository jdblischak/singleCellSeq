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
  -q --quantiles=<q>     Calculate the correlation for the genes separated by the provided
                         quantiles, e.g. -q .25 -q .75

Arguments:
  num_cells     number of single cells to subsample
  seed          seed for random number generator
  single        sample-by-gene matrix of single cell data
  bulk          gene-by-sample matrix of bulk cell data" -> doc

main <- function(num_cells, seed, single_fname, bulk_fname, individual = NULL,
                 good_cells = NULL, quantiles = NULL) {
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
  # Remove bulk samples
  single_cells <- single_cells[single_cells$well != "bulk", ]
  # Add rownames
  rownames(single_cells) <- paste(single_cells$individual, single_cells$batch,
                                  single_cells$well, sep = ".")
  # Remove meta-info cols
  single_cells <- single_cells[, grepl("ENSG", colnames(single_cells)) |
                                 grepl("ERCC", colnames(single_cells))]
  # Transpose
  single_cells <- t(single_cells)
  # Fix ERCC names
  rownames(single_cells) <- sub(pattern = "\\.", replacement = "-",
                                rownames(single_cells))
  # Keep only good quality cells
  if (!is.null(good_cells)) {
    assert("File with list of good quality cells exists.",
           file.exists(good_cells))
    good_cells_list <- scan(good_cells, what = "character", quiet = TRUE)
    good_cells_list <- substr(good_cells_list, start = 3, stop = 13)
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
  bulk_cells <- read.table(bulk_fname, row.names = 1, header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE)
  # Filter by individual
  if (!is.null(individual)) {
    bulk_cells <- bulk_cells[, grepl(individual, colnames(bulk_cells))]
  }
  # Remove single cells
  bulk_cells <- bulk_cells[, grepl("bulk", colnames(bulk_cells))]
  # Calculate cpm
  bulk_cells <- cpm(bulk_cells)

#   bulk_cells[1:10, 1:10]
#   dim(bulk_cells)
#   head(bulk_cells)

  # Filter genes
  # Remove unexpressed
  bulk_cells <- bulk_cells[rowSums(bulk_cells) > 0, ]
  # Remove the bottom 25%
  bulk_cells <- bulk_cells[rowMeans(bulk_cells) > quantile(rowMeans(bulk_cells), .25), ]
  # Filter genes in single cells
  single_cells <- single_cells[rownames(single_cells) %in% rownames(bulk_cells), ]

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
       quantiles = as.numeric(opts$quantiles))
} else if (interactive() & getOption('run.main', default = TRUE)) {
  # what to do if interactively testing
  main(num_cells = 20,
       seed = 1,
       single_fname = "/mnt/gluster/data/internal_supp/singleCellSeq/subsampled/molecule-counts-200000.txt",
       bulk_fname = "~/singleCellSeq/data/reads.txt",
       individual = "19098",
       good_cells = "../data/quality-single-cells.txt",
       quantiles = c(.25, .5, .75))
}
