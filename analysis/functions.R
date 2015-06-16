

run_pca <- function(x, retx = TRUE, center = TRUE, scale = TRUE) {
  pca <- prcomp(t(x), retx = TRUE, center = center, scale. = scale)
  variances <- pca$sdev^2
  explained <- variances / sum(variances)
  return(list(PCs = pca$x, explained = explained))
}


plot_pca <- function(x, pcx = 1, pcy = 2, explained = NULL, color = NULL,
                     shape = NULL, size = NULL) {
  library("ggplot2")
#   library("testit")
#   assert("PC and metadata  matrix have same number of rows.",
#          nrow(pc) == nrow(metadata))
  plot_data <- cbind(x, color, shape, size)
  plot_data <- as.data.frame(plot_data)
  p <- ggplot(plot_data, aes_string(x = paste0("PC", pcx),
                                    y = paste0("PC", pcy))) +
    geom_point()
  if (!is.null(color)) {
    p <- p + geom_point(aes(color = as.factor(color)))
  }
  if (!is.null(shape)) {
    p <- p + geom_point(aes(shape = as.factor(shape)))
  }


    # geom_point()
    # geom_point(aes(col = col, shape = shape, size = size))
  print(p)
}
