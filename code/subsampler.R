#!/usr/bin/env Rscript

# Run -h for command-line options.

library("testit")
suppressPackageStartupMessages(library("docopt"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("edgeR"))

"Calculate statistics on the subsampled data.

Usage:
subsampler.R [options] <num_cells> <seed> <single_sub> <single_full> <ercc> <bulk>

Options:
  -h --help              Show this screen.
  --individual=<ind>     Only use data from ind, e.g. NA19098
  --replicate=<rep>      Only use data from rep, e.g. r1
  --good_cells=<file>    A 1-column file with the names of good quality cells to maintain
  --keep_genes=<file>    A 1-column file with the names of genes to maintain
  -l --lower_q=<l>       The lower quantile cutoff to filter genes based on
                         expression level [default: 0]
  -u --upper_q=<u>       The upper quantile cutoff to filter genes based on
                         expression level [default: 1]
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
                 good_cells = NULL, keep_genes = NULL, lower_q = 0, upper_q = 1,
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

  # If an output filename is specified, save the PDF using the same name with
  # .pdf added as the file extension.
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

  # Split endogenous and ERCC genes
  ensg_index <- grepl("ENSG", rownames(bulk_cells))
  ercc_index <- grepl("ERCC", rownames(bulk_cells))
  assert("Both endogenous Ensembl and ERCC genes available for analysis",
         sum(ensg_index) > 0, sum(ercc_index) > 0)

  # Filter genes based on quantile cutoffs
  assert("Quantile cutoffs for gene expression level are valid input",
         is.numeric(lower_q), is.numeric(upper_q), lower_q >= 0, upper_q <= 1,
         upper_q > lower_q)
  results$lower_q <- lower_q
  results$upper_q <- upper_q
  # Use the rowSums of the counts in the bulk samples to order the genes by
  # expression level. This is sufficient because it maintains the relative gene
  # ordering.
  #cor(rowSums(bulk_cells), rowMeans(cpm(bulk_cells)), method = "spearman")
  relative_gene_exp <- rowSums(bulk_cells)
  # endogenous
  results$available_ensg <- length(relative_gene_exp[ensg_index])
  q_cutoffs_ensg <- quantile(relative_gene_exp[ensg_index],
                             probs = c(lower_q, upper_q))
  q_exp_filter_ensg <- relative_gene_exp >= q_cutoffs_ensg[1] &
                       relative_gene_exp <= q_cutoffs_ensg[2]
  results$used_ensg <- sum(ensg_index & q_exp_filter_ensg)
  # ercc
  results$available_ercc <- length(relative_gene_exp[ercc_index])
  q_cutoffs_ercc <- quantile(relative_gene_exp[ercc_index],
                             probs = c(lower_q, upper_q))
  q_exp_filter_ercc <- relative_gene_exp >= q_cutoffs_ercc[1] &
                       relative_gene_exp <= q_cutoffs_ercc[2]
  results$used_ercc <- sum(ercc_index & q_exp_filter_ercc)
  # Perform filtering
  single_cells_sub_ensg <- single_cells_sub[ensg_index & q_exp_filter_ensg, , drop = FALSE]
  single_cells_sub_ercc <- single_cells_sub[ercc_index & q_exp_filter_ercc, , drop = FALSE]
  single_cells_full_ensg <- single_cells_full[ensg_index & q_exp_filter_ensg, , drop = FALSE]
  single_cells_full_ercc <- single_cells_full[ercc_index & q_exp_filter_ercc, , drop = FALSE]
  bulk_cells_ensg <- bulk_cells[ensg_index & q_exp_filter_ensg, , drop = FALSE]
  bulk_cells_ercc <- bulk_cells[ercc_index & q_exp_filter_ercc, , drop = FALSE]
  assert("Correct number of genes after quantile filtering",
         nrow(single_cells_sub_ensg) == results$used_ensg,
         nrow(single_cells_sub_ercc) == results$used_ercc,
         nrow(single_cells_full_ensg) == results$used_ensg,
         nrow(single_cells_full_ercc) == results$used_ercc,
         nrow(bulk_cells_ensg) == results$used_ensg,
         nrow(bulk_cells_ercc) == results$used_ercc)

  # Some diagnostic plots exploring the total counts for endogenous and ERCC
  # genes between the subsampled and original full single cell data sets.
  if (diagnose) {
    totals_endo_single_sub <- colSums(single_cells_sub_ensg)
    totals_ercc_single_sub <- colSums(single_cells_sub_ercc)
    totals_endo_single_full <- colSums(single_cells_full_ensg)
    totals_ercc_single_full <- colSums(single_cells_full_ercc)
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
  subsample_index <- sample(1:ncol(single_cells_sub), size = num_cells)
  single_cells_sub_ensg <- single_cells_sub_ensg[, subsample_index]
  single_cells_sub_ercc <- single_cells_sub_ercc[, subsample_index]

  # Calculate mean gene expression:
  # For single cells, sum the counts across the single cells and then calculate
  # log2 cpm
  single_cells_cpm_mean_ensg <- calc_mean_single_cell(single_cells_sub_ensg)
  single_cells_cpm_mean_ercc <- calc_mean_single_cell(single_cells_sub_ercc)
  # For bulk samples, calculate cpm for each replicate and then calculate the
  # mean across the replicates
  bulk_cells_cpm_mean_ensg <- calc_mean_bulk_cell(bulk_cells_ensg)
  bulk_cells_cpm_mean_ercc <- calc_mean_bulk_cell(bulk_cells_ercc)

  # Calculate correlation between single cells and bulk samples
  results$pearson_ensg <- calc_cor(single_cells_cpm_mean_ensg, bulk_cells_cpm_mean_ensg,
                                   diagnose = diagnose, method = "pearson",
                                   prefix = "Correlation to bulk")
  results$pearson_ercc <- calc_cor(single_cells_cpm_mean_ercc, bulk_cells_cpm_mean_ercc,
                                   diagnose = diagnose, method = "pearson",
                                   prefix = "Correlation to bulk")
  results$spearman_ensg <- calc_cor(single_cells_cpm_mean_ensg, bulk_cells_cpm_mean_ensg,
                                    diagnose = FALSE, method = "spearman")
  results$spearman_ercc <- calc_cor(single_cells_cpm_mean_ercc, bulk_cells_cpm_mean_ercc,
                                    diagnose = FALSE, method = "spearman")

  # Detect number of expressed genes
  min_count = 1
  min_cells = 1
  detected_ensg_index <- apply(single_cells_sub_ensg, 1, detect_expression,
                               min_count = min_count, min_cells = min_cells)
  results$detected_ensg <- sum(detected_ensg_index)
  detected_ercc_index <- apply(single_cells_sub_ercc, 1, detect_expression,
                               min_count = min_count, min_cells = min_cells)
  results$detected_ercc <- sum(detected_ercc_index)

  # Detect number of reads/molecules
  # Caculate mean number of total counts, using only genes which meet the
  # criteria for detection.
  single_cells_sub_ensg_detected <- single_cells_sub_ensg[detected_ensg_index, , drop = FALSE]
  results$mean_counts_ensg <- mean(colSums(single_cells_sub_ensg_detected))
  single_cells_sub_ercc_detected <- single_cells_sub_ercc[detected_ercc_index, , drop = FALSE]
  results$mean_counts_ercc <- mean(colSums(single_cells_sub_ercc_detected))

  # Calculate variance
  var_full <- apply(single_cells_full_ensg, 1, var)
  var_full_log <- log(var_full + 0.25)
  var_sub <- apply(single_cells_sub_ensg, 1, var)
  var_sub_log <- log(var_sub + 0.25)

  results$var_pearson <- calc_cor(var_full_log, var_sub_log, method = "pearson",
                                  diagnose = diagnose, prefix = "variance")
  results$var_spearman <- calc_cor(var_full_log, var_sub_log, method = "spearman",
                                   diagnose = FALSE)

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
  x <- fread(fname, data.table = FALSE)
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
# prefix - a character prepended to the title of the diagnostic plot
#
# Returns a numeric vector of mean counts per million
#
calc_cor <- function(x, y, method = "pearson", diagnose = FALSE,
                          digits = 5, prefix = "") {
  assert("Proper input", length(x) == length(y), is.numeric(x), is.numeric(y),
         is.character(prefix))
  correlation <- cor(x, y, method = method)
  if (diagnose) {
    model <- lm(y ~ x)
    plot(x, y)
    abline(0, 1, col = "red")
    abline(model, col = "blue")
    title(main = sprintf("%s\tr: %.4f", prefix, correlation),
          sub = "Red: 1-1    Blue: Best fit")
  }
  assert("Proper output", correlation >= -1, correlation <= 1,
         length(correlation) == 1)
  correlation <- round(correlation, digits = digits)
  return(correlation)
}

# Detect which genes are expressed.
#
# x - vector of counts
# min_count - minumum required observations per cell
# min_cells - minimum required number of cells
#
# Returns TRUE if gene expressed, FALSE otherwise.
#
detect_expression <- function(x, min_count, min_cells) {
  assert("Proper input", is.numeric(x), is.numeric(min_count),
         is.numeric(min_cells), x >= 0, min_count > 0,
         length(min_count) == 1, length(min_cells) == 1)
  expressed <- sum(x >= min_count)
  return(expressed >= min_cells)
}

assert("Detect expression function works properly.",
       detect_expression(c(0, 0, 1, 1), 1, 2) == TRUE,
       detect_expression(c(0, 0, 0, 1), 1, 2) == FALSE,
       detect_expression(c(0, 0, 1, 1), 1, 3) == FALSE,
       detect_expression(c(0, 1, 1, 1), 1, 3) == TRUE,
       detect_expression(c(0, 2, 3, 1), 2, 2) == TRUE,
       detect_expression(c(0, 1, 3, 1), 2, 2) == FALSE)

# Input parameters -------------------------------------------------------------

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
       lower_q = as.numeric(opts$lower_q),
       upper_q = as.numeric(opts$upper_q),
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
       lower_q = 0.25,
       upper_q = 0.75,
       diagnose = TRUE)
}
