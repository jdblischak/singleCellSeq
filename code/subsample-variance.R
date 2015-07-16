#!/usr/bin/env Rscript

"Compare the variance estimates of subsamples to all cells.

Usage:
subsample-variance.R [options] <num_cells> <seed> <exp>

Options:
  -h --help              Show this screen.
  --individual=<ind>     Only use data from ind, e.g. 19098
  --min_count=<x>        The minimum count required for inclusion [default: 1]
  --min_cells=<x>        The minimum number of cells required for inclusion [default: 1]
  --good_cells=<file>    A 1-column file with the names of good quality cells to maintain

Arguments:
  num_cells     number of single cells to subsample
  seed          seed for random number generator
  exp           sample-by-gene expression matrix" -> doc

suppressMessages(library("docopt"))
library("testit")

main <- function(num_cells, seed, exp_fname, individual = NULL, min_count = 1,
                 min_cells = 1, good_cells = NULL) {
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
                                 grepl("ERCC", colnames(exp_dat)), drop = FALSE]
  # Transpose
  exp_dat <- t(exp_dat)
  # Fix ERCC names
  rownames(exp_dat) <- sub(pattern = "\\.", replacement = "-",
                                rownames(exp_dat))
  # Keep only good quality cells
  if (!is.null(good_cells)) {
    assert("File with list of good quality cells exists.",
           file.exists(good_cells))
    good_cells_list <- scan(good_cells, what = "character", quiet = TRUE)
    good_cells_list <- substr(good_cells_list, start = 3, stop = 13)
    exp_dat <- exp_dat[, colnames(exp_dat) %in% good_cells_list, drop = FALSE]
    assert("There are quality cells to perform the analysis.",
           ncol(exp_dat) > 0)
  }

  # Calculate total number of cells in full data set
  total_cells <- ncol(exp_dat)

  # Subsample number of single cells
  if (total_cells < num_cells) {
    cat(sprintf("%d\t%d\t%d\tNA\tNA\n", num_cells, seed, total_cells))
    return(invisible())
  }
  set.seed(seed)
  sub_indices <- sample(1:ncol(exp_dat), size = num_cells)

  # Filter to only include expressed genes
  expressed <- apply(exp_dat, 1, detect_expression, min_count = min_count,
                     min_cells = min_cells)
  exp_dat <- exp_dat[expressed, ]

  # Calculate variance
  var_full <- apply(exp_dat, 1, var)
  var_full <- log(var_full + 0.25)
  var_sub <- apply(exp_dat[, sub_indices, drop = FALSE], 1, var)
  var_sub <- log(var_sub + 0.25)

  # Calculate Pearson correlation coefficient
  r <- cor(var_full, var_sub)

  # Calculate root-mean-square error (RMSE)
  rmse <- calc_rmse(var_full, var_sub)

  # Output
  cat(sprintf("%d\t%d\t%d\t%.2f\t%.2f\n", num_cells, seed, total_cells, r, rmse))
}

detect_expression <- function(x, min_count, min_cells) {
  # x - vector of counts
  # min_count - minumum required observations per cell
  # min_cells - minimum required number of cells
  expressed <- sum(x >= min_count)
  return(expressed >= min_cells)
}

library("testit")
assert("Detect expression function works properly.",
       detect_expression(c(0, 0, 1, 1), 1, 2) == TRUE,
       detect_expression(c(0, 0, 0, 1), 1, 2) == FALSE,
       detect_expression(c(0, 0, 1, 1), 1, 3) == FALSE,
       detect_expression(c(0, 1, 1, 1), 1, 3) == TRUE,
       detect_expression(c(0, 2, 3, 1), 2, 2) == TRUE,
       detect_expression(c(0, 1, 3, 1), 2, 2) == FALSE)

calc_rmse <- function(x, y) {
  # Calculate the root-mean-square error of the two vectors
  assert("Vectors of same length",
         length(x) == length(y))
  sse <- (x - y)^2
  mse <- mean(sse)
  rmse <- sqrt(mse)
  return(rmse)
}

assert("RMSE calculated correctly",
       calc_rmse(c(1, 0, -1), c(0, -1, 1)) == sqrt(2))

if (!interactive() & getOption('run.main', default = TRUE)) {
  opts <- docopt(doc)
  # str(opts)
  main(num_cells = as.numeric(opts$num_cells),
       seed = as.numeric(opts$seed),
       exp_fname = opts$exp,
       min_count = as.numeric(opts$min_count),
       min_cells = as.numeric(opts$min_cells),
       individual = opts$individual,
       good_cells = opts$good_cells)
} else if (interactive() & getOption('run.main', default = TRUE)) {
  # what to do if interactively testing
  main(num_cells = 20,
       seed = 1,
       exp_fname = "/mnt/gluster/data/internal_supp/singleCellSeq/subsampled/molecule-counts-200000.txt",
       individual = "19098",
       min_count = 10,
       min_cells = 5,
       good_cells = "../data/quality-single-cells.txt")
}
