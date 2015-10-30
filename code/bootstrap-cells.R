#' Bootstrap cell labels within each individual
#'
#' @param log2counts log2counts matrix of gene by cells
#' @param grouping_vector the grouping vector corresponds to variable of interest
#' @param num_bootstrap number of bootstrap samples
#'
#' @export
#'
#' @examples
#' bootstrap_cells()
#'
bootstrap_cells <- function(log2counts, grouping_vector, num_bootstrap) {
#   log2counts <- molecules_ENSG
#   grouping_vector <- anno_qc_filter$individual
#   num_permute <- 10
  num_cells <- dim(log2counts)[2]

  bootstrap_log2counts <- lapply(1:num_bootstrap, function(ii_group) {
    # creata a sequence of random numbers
    per_group <- lapply(1:3, function(ii_group) {
      ind_log2counts <- log2counts[ , grouping_vector == unique(grouping_vector)[ii_group]]
      num_cells <- ncol(ind_log2counts)
      bootstrap_data <- ind_log2counts[ , sample(1:num_cells)]
      bootstrap_data
    })
    per_group <- do.call(cbind, per_group)
    per_group
  })
  return(bootstrap_log2counts)
}




