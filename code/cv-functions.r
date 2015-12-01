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
                           subsample_size,
                           number_subsample)
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
                                          subsample_size + i)
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


#' Compute coefficient of variation 
#' 
#' Compute CV across cells belong to each level of the grouping variable.
#'
#' @param log2counts log2 count matrix of gene by cell.
#' @param grouping_vector Compute per gene CV for cells belonging to each level of
#'        the grouping vector.
#'
#' @export
#'
#' @examples
#' compute_cv()
#'
compute_cv <- function(log2counts, grouping_vector) {

  group_cv <- lapply( unique(grouping_vector), function(per_group) {
    # Convert log2cpm to counts
    counts_per_group <- 2^log2counts[ , grouping_vector == per_group ]
    mean_per_gene <- apply(counts_per_group, 1, mean, na.rm = TRUE)
    sd_per_gene <- apply(counts_per_group, 1, sd, na.rm = TRUE)
    cv_per_gene <- data.frame(mean = mean_per_gene,
                              sd = sd_per_gene,
                              cv = sd_per_gene/mean_per_gene,
                              group = rep(per_group, dim(counts_per_group)[1]) )
    rownames(cv_per_gene) <- rownames(counts_per_group)

    return(cv_per_gene)
  })
  names(group_cv) <- unique(grouping_vector)
  group_cv
}


#' Normalize coefficients of variation
#'
#' @param group_cv CVs per batch computed use compute_cv().
#' @param log2counts log2 count matrix of gene by cell.
#' 
#' @export
#'
#' @examples
#' normalize_cv()
#'
normalize_cv <- function(group_cv, log2counts, anno) {
  library(zoo)
  # Compute a data-wide coefficient of variation on counts.
  data_cv <- apply(2^log2counts, 1, sd)/apply(2^log2counts, 1, mean)

  # Order genes by mean expression levels
  order_gene <- order(apply(2^log2counts, 1, mean))

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
  names(roll_medians) <- rownames(log2counts)[order_gene]

  # Order rolling medians according to the count matrix
  reorder_gene <- match(rownames(log2counts), names(roll_medians) )
  roll_medians <- roll_medians[ reorder_gene ]
  stopifnot( all.equal(names(roll_medians), rownames(log2counts) ) )

  group_cv_adj <- lapply(1:length(group_cv), function(ii_batch) {
    # Adjusted coefficient of variation on log10 scale
    log10cv2_adj <- log10(group_cv[[ii_batch]]$cv^2) - roll_medians
    # combine the adjusted cv with the unadjusted cv
    data.frame(group_cv[[ii_batch]],
               log10cv2_adj = log10cv2_adj )
  })
  names(group_cv_adj) <- names(group_cv)
  group_cv_adj
}
