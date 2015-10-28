#' Permute cell labels
#'
#' @param log2counts log2counts matrix of gene by cells
#' @param grouping_vector the grouping vector corresponds to variable of interest
#' @param num_permute number of permuted samples
#'
#' @export
#'
#' @examples
#' permuate_cells()
permute_cells <- function(log2counts, grouping_vector, num_permute) {
#   log2counts <- molecules_ENSG
#   grouping_vector <- anno_qc_filter$individual
#   num_permute <- 10
  num_cells <- dim(log2counts)[2]

  permuted_log2counts <- lapply(1:num_permute, function(ii_permute) {
    # creata a sequence of random numbers
    perm_labels <- sample(num_cells, replace = TRUE)

    # reorder columns (cells) of the data matrix
    # according to perm_labels
    perm_data <- log2counts[ , perm_labels]

    # label the columns using the original individual labels
    # now the cells for individual A in the permuted data
    # can come from indivdual A, B, or C
    colnames(perm_data) <- grouping_vector
    perm_data
  })
  return(permuted_log2counts)
}




