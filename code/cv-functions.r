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

  log2counts <- as.matrix(log2counts)

  groups <- c(unique(grouping_vector), "all")

  group_cv <- lapply( groups, function(per_group) {
    if(per_group != "all") {
      # Convert log2cpm to counts
      counts_per_group <- 2^log2counts[ , grouping_vector == per_group ]
    }
    if (per_group == "all") {
      # Convert log2cpm to counts
      counts_per_group <- 2^log2counts
    }
    mean_per_gene <- apply(counts_per_group, 1, mean, na.rm = TRUE)
    sd_per_gene <- apply(counts_per_group, 1, sd, na.rm = TRUE)
    cv_per_gene <- data.frame(mean = mean_per_gene,
                              sd = sd_per_gene,
                              cv = sd_per_gene/mean_per_gene,
                              group = rep(per_group, dim(counts_per_group)[1]) )
    rownames(cv_per_gene) <- rownames(counts_per_group)

    return(cv_per_gene)
  })
  names(group_cv) <- groups
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




#' Compute CV of expressed cells
#'
#' @param molecules_filter molecule count after collision correction.
#' @param molecules_final molecule count after standardization and batch correction.
#'
#' @export
#'
#' @examples
#' expressed_cv()
#'
#' #aa <- expressed_cv(molecules_filter,
#' #                   molecules_final,
#' #                   anno_filter$individual)
#'
compute_expressed_cv <- function(molecules_filter,
                         molecules_final,
                         grouping_vector) {
  groups <- c(unique(grouping_vector), "all")

  # Set all zero count cells to be NA
  molecules_filter_df <- molecules_filter[match(rownames(molecules_final),
                                             rownames(molecules_filter)), ]
  molecules_filter_df[which(molecules_filter_df == 0, arr.ind = TRUE)] <- NA
  molecules_filter_df[which(molecules_filter_df != 0, arr.ind = TRUE)] <- 1

  molecules_final_df <- molecules_final

  expressed_cv <- lapply( groups, function(ind) {
    if(ind != "all") {
      temp_final_ind <- molecules_final_df[, anno_filter$individual == ind]
      temp_filter_ind <- molecules_filter_df[, anno_filter$individual == ind]
    }
    if(ind == "all") {
      temp_final_ind <- molecules_final_df
      temp_filter_ind <- molecules_filter_df
    }
  # exclude genes with zero counts in all cells
    which_include <- which(rowSums(is.na(temp_filter_ind)) != NCOL(temp_filter_ind))
    temp_filter_ind <- temp_filter_ind[which_include, ]
    temp_final_ind <- temp_final_ind[which_include, ]

    temp_ind <- 2^(temp_final_ind*temp_filter_ind)
    gene_expr <- data.frame(
      expr_mean = rowMeans(as.matrix(temp_ind), na.rm = TRUE),
      expr_var = matrixStats::rowVars(as.matrix(temp_ind), na.rm = TRUE),
      expr_cv = sqrt(matrixStats::rowVars(as.matrix(temp_ind), na.rm = TRUE))/rowMeans(as.matrix(temp_ind), na.rm = TRUE)
    )

    rownames(gene_expr) <- rownames(temp_final_ind)
    return(gene_expr)
  })
  names(expressed_cv) <- groups

  # filter genes that are present in all three individuals
  which_genes_all <- Reduce("intersect",
                             lapply(expressed_cv, rownames))

  expressed_cv_filter <- lapply(expressed_cv, function(x) {
      x[which(rownames(x) %in% which_genes_all), ]
  })
  names(expressed_cv_filter) <- names(expressed_cv)


  return(expressed_cv_filter)
}





#' Normalize coefficients of variation with summary statistics
#'
#' @param expressed_cv a list. Each list is a data.frame of
#'        gene mean expression, gene expression cv, and
#'        gene expression variance. The last list is a data.frame of
#'        the same summary statistics across all individuals.
#'
#' @export
#'
#' @examples
#' normalize_cv_input()
#'
normalize_cv_input <- function(expressed_cv,
                               grouping_vector) {
  groups <- unique(grouping_vector)
  ## normalize CV
  library(zoo)
  # order genes by overall mean expression level
  expressed_gene_order <- order(expressed_cv[["all"]]$expr_mean)

  # Rolling medians of log10 squared CV by mean expression levels
  roll_medians <- suppressWarnings(
    rollapply( log10(expressed_cv[["all"]]$expr_cv^2)[expressed_gene_order],
               width = 50, by = 25,
               FUN = median, fill = list("extend", "extend", "NA") )
  )
  ii_na <- which( is.na(roll_medians) )
  roll_medians[ii_na] <- median(
    log10(expressed_cv[["all"]]$expr_cv^2)[expressed_gene_order][ii_na] )
  names(roll_medians) <- rownames(expressed_cv[[1]])[expressed_gene_order]

  # Order rolling medians according to the count matrix
  reorder_gene <- match(rownames(expressed_cv[[1]]), names(roll_medians) )
  roll_medians <- roll_medians[ reorder_gene ]
  stopifnot( all.equal(names(roll_medians), rownames(expressed_cv[[1]]) ) )

  expressed_dm <- do.call(cbind,
                          lapply(groups, function(ind) {
                            # Adjusted coefficient of variation on log10 scale
                            dm <- log10(expressed_cv[[ind]]$expr_cv^2) - roll_medians
                            # combine the adjusted cv with the unadjusted cv
                            return(dm)
                          }) )
  colnames(expressed_dm) <- groups
  expressed_dm <- data.frame(expressed_dm)
  return(expressed_dm)
}




#' Normalize coefficients of variation with summary statistics
#'
#' @param expressed_cv a list. Each list is a data.frame of
#'        gene mean expression, gene expression cv, and
#'        gene expression variance. The last list is a data.frame of
#'        the same summary statistics across all individuals.
#'
#' @export
#'
#' @examples
#' normalize_cv_input()
#'
normalize_cv_input <- function(expressed_cv,
                               grouping_vector) {
  groups <- unique(grouping_vector)
  ## normalize CV
  library(zoo)
  # order genes by overall mean expression level
  expressed_gene_order <- order(expressed_cv[["all"]]$expr_mean)

  # Rolling medians of log10 squared CV by mean expression levels
  roll_medians <- suppressWarnings(
    rollapply( log10(expressed_cv[["all"]]$expr_cv^2)[expressed_gene_order],
               width = 50, by = 25,
               FUN = median, fill = list("extend", "extend", "NA") )
  )
  ii_na <- which( is.na(roll_medians) )
  roll_medians[ii_na] <- median(
    log10(expressed_cv[["all"]]$expr_cv^2)[expressed_gene_order][ii_na] )
  names(roll_medians) <- rownames(expressed_cv[[1]])[expressed_gene_order]

  # Order rolling medians according to the count matrix
  reorder_gene <- match(rownames(expressed_cv[[1]]), names(roll_medians) )
  roll_medians <- roll_medians[ reorder_gene ]
  stopifnot( all.equal(names(roll_medians), rownames(expressed_cv[[1]]) ) )

  expressed_dm <- do.call(cbind,
                          lapply(groups, function(ind) {
                            # Adjusted coefficient of variation on log10 scale
                            dm <- log10(expressed_cv[[ind]]$expr_cv^2) - roll_medians
                            # combine the adjusted cv with the unadjusted cv
                            return(dm)
                          }) )
  colnames(expressed_dm) <- groups
  expressed_dm <- data.frame(expressed_dm,
                             all = log10(expressed_cv[["all"]]$expr_cv^2) - roll_medians)

  return(expressed_dm)
}
