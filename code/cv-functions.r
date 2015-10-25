####----- Coefficient of variation related functions -----####

#' Subsample molecule count (or read count) data
#'
#' @param counts Count matrix of gene by cell.
#' @param anno Annotation matrix labeling each cell according to individual,
#'        batch, replicate, well ID, etc.
#' @param subsample_size Number of cells in each subsample.
#' @param number_subsample Number of subsamples to be taken from each individual.
#'
#' @export
#'
#' @examples
#' make_subsample()

make_subsample <- function(counts, anno,
                           subsample_size = 20,
                           number_subsample = 5)
{
  ind_subsample <- lapply( unique(anno$individual),
                           function(per_individual) {
                             # Indicator variable for the selected individual
                             ii_individual <- anno$individual == per_individual
                             # Total number of cells for the selected individual
                             nn_individual <- sum(ii_individual)
                             # Annotation matrix for the selected individual
                             anno_individual <- anno[ii_individual, ]
                             # Normalized counts for the selected individual
                             counts_individual <- counts[ , ii_individual]
                             # Subsample counts
                             counts_subsample <- lapply(1:number_subsample, function(i) {
                               set.seed(which(unique(anno$individual) == per_individual)*1000 +
                                          subsample_size + number_subsample)
                               ii_sample <- sample(nn_individual, subsample_size, replace = FALSE)
                               counts_subsample <- counts_individual[ , ii_sample]
                               colnames(counts_subsample) <- rep(paste0("r.", i), ncol(counts_subsample) )
                               counts_subsample
                             })
                             counts_subsample <- do.call(cbind, counts_subsample)

                             # Subsample annotation
                             anno_subsample <- data.frame(individual = rep(per_individual, number_subsample*subsample_size),
                                                          replicate = colnames(counts_subsample) )
                             anno_subsample$batch <- with(anno_subsample, paste0(individual, replicate))

                             list(counts_subsample = counts_subsample,
                                  anno_subsample = anno_subsample)
                           })

  counts_subsample <- do.call(cbind, lapply(ind_subsample, "[[", 1))
  anno_subsample <- do.call(rbind, lapply(ind_subsample, "[[", 2))

  list(counts_subsample = counts_subsample,
       anno_subsample = anno_subsample)
}


#' Compute coefficient of variation for each batch (per individual, per replicate)
#'
#' @param counts Count matrix of gene by cell.
#' @param anno Annotation matrix labeling each cell according to individual,
#'        batch, replicate, well ID, etc.
#'
#' @export
#'
#' @examples
#' compute_cv()
#'
compute_cv <- function(counts, anno) {
  # Compute CV for each batch
  batch_cv <- lapply( unique(anno$batch), function(per_batch) {
    # Convert log2cpm to counts
    counts_per_batch <- 2^counts[ , anno$batch == per_batch ]
    anno_per_batch <- anno[ anno$batch == per_batch, ]
    mean_per_gene <- apply(counts_per_batch, 1, mean, na.rm = TRUE)
    sd_per_gene <- apply(counts_per_batch, 1, sd, na.rm = TRUE)
    cv_per_gene <- data.frame(mean = mean_per_gene,
                              sd = sd_per_gene,
                              cv = sd_per_gene/mean_per_gene,
                              individual = unique(anno_per_batch$individual),
                              replicate = unique(anno_per_batch$replicate),
                              batch = unique(anno_per_batch$batch) )
    rownames(cv_per_gene) <- rownames(counts_per_batch)

    return(cv_per_gene)
  })
  names(batch_cv) <- unique(anno$batch)
  batch_cv
}

#' Subsample molecule count (or read count) data
#'
#' @param batch_cv CVs per batch computed use compute_cv().
#' @param counts Count matrix of gene by cell.
#' @param anno Annotation matrix labeling each cell according to individual,
#'        batch, replicate, well ID, etc.
#'
#' @export
#'
#' @examples
#' normalize_cv()
#'
normalize_cv <- function(batch_cv, counts, anno) {
  library(zoo)
  # Compute a data-wide coefficient of variation on counts.
  data_cv <- apply(2^counts, 1, sd)/apply(2^counts, 1, mean)

  # Order genes by mean expression levels
  order_gene <- order(apply(2^counts, 1, mean))

  # Rolling medians of log10 squared CV by mean expression levels
  # Avoid warning message introduced by NA in the rollapply results,
  # which are introduced by coercion
  roll_medians <- suppressWarnings(
    rollapply( log10(data_cv^2)[order_gene],
               width = 50, by = 25,
               FUN = median, fill = list("extend", "extend", "NA") )
  )
  ii_na <- which( is.na(roll_medians) )
  roll_medians[ii_na] <- median( log10(data_cv^2)[order_gene][ii_na] )
  names(roll_medians) <- rownames(molecules_ENSG)[order_gene]

  # Order rolling medians according to the count matrix
  reorder_gene <- match(rownames(counts), names(roll_medians) )
  roll_medians <- roll_medians[ reorder_gene ]
  stopifnot( all.equal(names(roll_medians), rownames(counts) ) )

  batch_cv_adj <- lapply(1:length(batch_cv), function(ii_batch) {
    # Adjusted coefficient of variation on log10 scale
    log10cv2_adj <- log10(batch_cv[[ii_batch]]$cv^2) - roll_medians
    # combine the adjusted cv with the unadjusted cv
    data.frame(batch_cv[[ii_batch]],
               log10cv2_adj = log10cv2_adj )
  })
  names(batch_cv_adj) <- names(batch_cv)
  batch_cv_adj
}
