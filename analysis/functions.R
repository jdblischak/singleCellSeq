

run_pca <- function(x, retx = TRUE, center = TRUE, scale = TRUE) {
  pca <- prcomp(t(x), retx = TRUE, center = center, scale. = scale)
  variances <- pca$sdev^2
  explained <- variances / sum(variances)
  return(list(PCs = pca$x, explained = explained))
}


# Plot PCA results.
#
# Args:
#   x: numeric matrix of PCs
#   pcx: PC to plot on x-axis (default: 1)
#   pcy: PC to plot on y-axis (default: 2)
#   explained: numeric vector of fractions of variance explained by each PC
#   metadata: data frame or matrix that contains the metadata used to annotate
#             the plot
#   color, shape, size: column name of metadata used to pass column to ggplot
#                       aesthetic
#   factors: character vector which contains the column names of metadata that
#            need to be explcitly converted to a factor
plot_pca <- function(x, pcx = 1, pcy = 2, explained = NULL, metadata = NULL,
                     color = NULL, shape = NULL, size = NULL, factors = NULL) {
  library("ggplot2")
  library("testit")

  # Prepare data frame to pass to ggplot
  if (!is.null(metadata)) {
    assert("PC and metadata have same number of rows.",
           nrow(x) == nrow(metadata))
    plot_data <- cbind(x, metadata)
    plot_data <- as.data.frame(plot_data)
    # Convert numeric factors to class "factor"
    for (f in factors) {
      plot_data[, f] <- as.factor(plot_data[, f])
    }
  } else {
    plot_data <- as.data.frame(x)
  }
  # Prepare axis labels
  if (!is.null(explained)) {
    assert("Number of PCs differs between x and explained.",
           length(explained) == ncol(x))
    xaxis <- sprintf("PC%d (%.2f%%)", pcx, round(explained[pcx] * 100, 2))
    yaxis <- sprintf("PC%d (%.2f%%)", pcy, round(explained[pcy] * 100, 2))
  } else {
     xaxis <- paste0("PC", pcx)
     yaxis <- paste0("PC", pcy)
  }
  # Plot
  p <- ggplot(plot_data, aes_string(x = paste0("PC", pcx),
                                    y = paste0("PC", pcy))) +
    geom_point(aes_string(color = color, shape = shape, size = size)) +
    labs(x = xaxis, y = yaxis)
  p
}
