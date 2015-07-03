#!/usr/bin/env Rscript

"Detect the number of expressed genes.

Usage:
detect-genes.R [options] <num_cells> <seed> <exp>

Options:
  -h --help              Show this screen.
  --individual=<ind>     Only use data from ind, e.g. 19098
  --count=<x>            The minimum count required for detection [default: 1]
  --fraction=<x>         The minimum fraction of cells required for detection [default: 0.5]

Arguments:
  num_cells     number of single cells to subsample
  seed          seed for random number generator
  exp           sample-by-gene expression matrix" -> doc

suppressMessages(library("docopt"))

main <- function(num_cells, seed, exp_fname, individual = NULL, count = 1,
                 fraction = 0.5) {
  # Load expression data
  exp_dat <- read.table(exp_fname, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)
  # Filter by individual
  if (!is.null(individual)) {
    exp_dat <- exp_dat[exp_dat$individual == individual, ]
  }
  # Remove bulk samples
  exp_dat <- exp_dat[exp_dat$well != "bulk", ]
  # Add rownames
  rownames(exp_dat) <- paste(exp_dat$individual, exp_dat$batch,
                                  exp_dat$well, sep = ".")
  # Remove meta-info cols
  exp_dat <- exp_dat[, grepl("ENSG", colnames(exp_dat)) |
                                 grepl("ERCC", colnames(exp_dat))]
  # Transpose
  exp_dat <- t(exp_dat)
  # Fix ERCC names
  rownames(exp_dat) <- sub(pattern = "\\.", replacement = "-",
                                rownames(exp_dat))

  # Subsample number of single cells
  if (ncol(exp_dat) < num_cells) {
    cat(sprintf("%d\t%d\tNA\n", num_cells, seed))
    quit()
  }
  set.seed(seed)
  exp_dat <- exp_dat[, sample(1:ncol(exp_dat), size = num_cells)]
#   exp_dat[1:10, 1:10]
#   dim(exp_dat)
  detected <- apply(exp_dat, 1, detect_expression, count = count, fraction = fraction)
  num_detected <- sum(detected)

  # Output
  cat(sprintf("%d\t%d\t%d\n", num_cells, seed, num_detected))
}

detect_expression <- function(x, count, fraction) {
  # x - vector of counts
  # count - minumum required observations per cell
  # fraction - minimum required fraction of cells
  expressed <- sum(x >= count)
  fraction_expressed <- expressed / length(x)
  return(fraction_expressed >= fraction)
}

library("testit")
assert("Detect expression function works properly.",
       detect_expression(c(0, 0, 1, 1), 1, 0.5) == TRUE,
       detect_expression(c(0, 0, 0, 1), 1, 0.5) == FALSE,
       detect_expression(c(0, 0, 1, 1), 1, 0.75) == FALSE,
       detect_expression(c(0, 1, 1, 1), 1, 0.75) == TRUE,
       detect_expression(c(0, 2, 3, 1), 2, 0.5) == TRUE,
       detect_expression(c(0, 1, 3, 1), 2, 0.5) == FALSE)

if (!interactive() & getOption('run.main', default = TRUE)) {
  opts <- docopt(doc)
  # str(opts)
  main(num_cells = as.numeric(opts$num_cells),
       seed = as.numeric(opts$seed),
       exp_fname = opts$exp,
       count <- as.numeric(opts$count),
       fraction <- as.numeric(opts$fraction),
       individual = opts$individual)
} else if (interactive() & getOption('run.main', default = TRUE)) {
  # what to do if interactively testing
  main(num_cells = 20,
       seed = 1,
       exp_fname = "/mnt/gluster/data/internal_supp/singleCellSeq/subsampled/molecule-counts-200000.txt",
       individual = "19098",
       count = 10,
       fraction = 0.5)
}
