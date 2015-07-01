#!/usr/bin/env Rscript

# Run -h for command-line options.

suppressMessages(library("docopt"))

"Correlate the expression in single cells to the bulk sample.

Usage:
correlate-single-cell-to-bulk.R [--individual=<ind>] <single> <bulk>

Options:
  -h --help              Show this screen.
  --individual=<ind>     Show version.

Arguments:
  single        sample-by-gene matrix of single cell data
  bulk          gene-by-sample matrix of bulk cell data" -> doc

main <- function(single_fname, bulk_fname, individual = NULL) {

  id <- "single-to-bulk-correlation"

  # Load single cell data
  single_cells <- read.table(single_fname,
                             header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # Filter by individual
  ind <- "19098"
  single_cells <- single_cells[single_cells$individual == ind, ]
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
  num_cells <- 20
  seed <- 1
  set.seed(seed)
  single_cells <- single_cells[, sample(1:ncol(single_cells), size = num_cells)]
  # Calculate cpm
  suppressPackageStartupMessages(library("edgeR"))
  single_cells <- cpm(single_cells)
#   single_cells[1:10, 1:10]
#   dim(single_cells)

  # Load bulk cell data
  bulk_cells <- read.table(bulk_fname, row.names = 1,
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # Filter by individual
  bulk_cells <- bulk_cells[, grepl(ind, colnames(bulk_cells))]
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
  library("testit")
  assert("Same number of genes in bulk and single cells.",
         nrow(bulk_cells) == nrow(single_cells))

  # Correlate
  r <- cor(rowMeans(single_cells), rowMeans(bulk_cells))

  # Output
  cat(sprintf("%s\t%d\t%d\t%f\n", ind, num_cells, seed, r),
      file = sprintf("%s-%d-%d.txt", id, num_cells, seed))

}


if (!interactive() & getOption('run.main', default = TRUE)) {
  opts <- docopt(doc)
  main(single_fname = opts$single,
       bulk_fname = opts$bulk)
} else if (interactive() & getOption('run.main', default = TRUE)) {
  # what to do if interactively testing
  main(single_fname = "/mnt/gluster/data/internal_supp/singleCellSeq/subsampled/molecule-counts-200000.txt",
       bulk_fname = "~/singleCellSeq/data/reads.txt")
}
