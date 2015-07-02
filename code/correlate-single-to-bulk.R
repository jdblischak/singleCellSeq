#!/usr/bin/env Rscript

# Run -h for command-line options.

suppressMessages(library("docopt"))

"Correlate the expression in single cells to the bulk sample.

Usage:
correlate-single-cell-to-bulk.R [--individual=<ind>] <num_cells> <seed> <single> <bulk>

Options:
  -h --help              Show this screen.
  --individual=<ind>     Only use data from ind, e.g. 19098

Arguments:
  num_cells     number of single cells to subsample
  seed          seed for random number generator
  single        sample-by-gene matrix of single cell data
  bulk          gene-by-sample matrix of bulk cell data" -> doc

main <- function(num_cells, seed, single_fname, bulk_fname, individual = NULL) {
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

  # Subsample number of single cells
  if (ncol(single_cells) < num_cells) {
    cat(sprintf("%d\t%d\tNA\n", num_cells, seed))
    quit()
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

  # Correlate
  r <- cor(rowMeans(single_cells), rowMeans(bulk_cells))

  # Output
  cat(sprintf("%d\t%d\t%f\n", num_cells, seed, r))
}


if (!interactive() & getOption('run.main', default = TRUE)) {
  opts <- docopt(doc)
  main(num_cells = as.numeric(opts$num_cells),
       seed = as.numeric(opts$seed),
       single_fname = opts$single,
       bulk_fname = opts$bulk,
       individual = opts$individual)
} else if (interactive() & getOption('run.main', default = TRUE)) {
  # what to do if interactively testing
  main(num_cells = 20,
       seed = 1,
       single_fname = "/mnt/gluster/data/internal_supp/singleCellSeq/subsampled/molecule-counts-200000.txt",
       bulk_fname = "~/singleCellSeq/data/reads.txt",
       individual = "19098")
}
