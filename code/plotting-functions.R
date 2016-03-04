####<----- Some plotting tools ----->####

#' Overlay density plots of log-transformed counts (for one gene)
#'
#' @param molecules Matrix of log-transformed counts. Sample by cell.
#' @param annotation Annotation matrix labeling each cell according to individual,
#'        batch, replicate, well ID, etc.
#' @param individual_label Individual (cell line) that we want to plot.
#' @param batches TRUE if we want to plot by batches.
#' @param which_gene Ensembl ID of the gene.
#' @param labels Subtitle of the plot.
#' @param gene_symbols Matrix listing each gene's Ensembl ID and HGNC symbol.
#'
#' @export
#'
#' @examples
#' make_subsample()
plot_density_overlay <- function(molecules, annotation,
                                individual_label = NULL, batches = NULL,
                                which_gene, labels,
                                xlims = NULL, ylims = NULL, gene_symbols) {
  if_present <- which(rownames(molecules) == which_gene)
  if(length(if_present) == 0) {
    stop("Gene not present in the data")
  }

  library(scales)
  library(broman)
  crayon <- brocolors("crayon")
  colors <- c("Sunset Orange", "Tropical Rain Forest", "Denim")

  if (!is.null(individual_label)) {
    annotation_df <- annotation[
      which(annotation$individual == individual_label), ]
    molecules_df <- molecules[ ,
                               which(annotation$individual == individual_label)]
    colors_df <- colors[which(unique(annotation$individual) == individual_label)]
    individuals <- individual_label
  } else {
    annotation_df <- annotation
    molecules_df <- molecules
    colors_df <- colors
    individuals <- unique(annotation_df$individual)
  }
  # compute density estimates
  dens <- lapply(individuals, function(ind) {
    df <- unlist(molecules_df[ rownames(molecules_df) == which_gene,
                               annotation_df$individual == ind])
    dens_try <- try(density(df, na.rm = TRUE), silent = TRUE)
    if(class(dens_try) == "try-error") {
      return(NULL)
    } else {
      return(dens_try)
    }
  })

  # find ranges of x and y-axis
  if (is.null(xlims)) xlims <- range(sapply(dens, function(obj) obj$x))
  if (is.null(ylims)) ylims <- range(sapply(dens, function(obj) obj$y))

  plot(dens[[1]],
       xlab = "log2 gene expression", main = "",
       ylab = "Density", axes = F, lwd = 0, xlim = xlims, ylim = ylims)
  for (i in 1:length(individuals)) {
    polygon(dens[[i]],
            col = alpha(crayon[colors_df[i]], .4),
            border = "grey40")
  }
  axis(1); axis(2)
  mtext(text = labels, side = 3)
  title(main = with(gene_symbols,
                    external_gene_name[which(ensembl_gene_id == which_gene)]) )

  #   # TBD
  #   if (!is.null(batches)) {
  #
  #     individuals <- unique(annotation$individual)
  #     dens <- lapply(1:length(individuals), function(per_individual) {
  #       which_individual <- annotation$individual == individuals[per_individual]
  #       annotation_sub <- annotation[which_individual, ]
  #       molecules_sub <- molecules_ENSG[ , which_individual]
  #       replicates <- unique(annotation_sub$replicate)
  #       dens_batch <- lapply(1:length(replicates), function(per_replicate) {
  #         which_replicate <- annotation_sub$replicate == replicates[per_replicate]
  #         density(unlist( molecules_sub[ rownames(molecules_ENSG) == which_gene,
  #                                        which_replicate] ) )
  #       })
  #     })
  #
  #     if (is.null(xlims)) {
  #       xlims <- range( c( sapply(dens, function(obj_individual) {
  #         c( sapply(obj_individual, function(obj) {
  #           range(obj$x)
  #         }) )
  #       }) ) )
  #     }
  #     if (is.null(ylims)) {
  #       ylims <- range( c( sapply(dens, function(obj_individual) {
  #         c( sapply(obj_individual, function(obj) {
  #           range(obj$y)
  #         }) )
  #       }) ) )
  #     }
  #
  #     colors <- c("Sunset Orange", "Tropical Rain Forest", "Denim")
  #     for (i in 1:length(dens)) {
  #
  #       if (i == 1) col <- crayon["Sunset Orange"]
  #       if (i == 2) col <- crayon["Tropical Rain Forest"]
  #       if (i == 3) col <- crayon["Denim"]
  #
  #       plot(dens[[i]][[1]],
  #            xlab = "log2 gene expression", main = "",
  #            ylab = "Density", axes = F, lwd = 0, xlim = xlims, ylim = ylims)
  #       for (j in 1:length(dens[[i]])) {
  #         polygon(dens[[i]][[j]],
  #                 col = alpha(col, .4),
  #                 border = "grey40")
  #       }
  #     }

  axis(1); axis(2)
  mtext(text = labels, side = 3)
  title(main = with(gene_symbols,
                    external_gene_name[which(ensembl_gene_id == which_gene)]) )
}




#' Coefficient of variation versus mean plot for the expressed cells
#'
#' @param expr_mean Gene mean expression vector
#' @param exprs_cv GEne coefficient of variation vector
#' @param main Plot label
#' @param ylab y-axis label
#'
#' @export
#'
#' @examples
#' plot_poisson_cv_expressed()

plot_poisson_cv_expressed <-
  function(expr_mean = NULL,
           exprs_cv = NULL,
           main, ylab = NULL) {

    # define colors
    cbPalette <- c("#999999", "#0000FF", "#56B4E9", "#009E73",
                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

    # exclude genes with missing mean or variance
    gene_valid <- intersect(which(!is.na(expr_mean)), which(!is.na(exprs_cv)) )

    # defnine poisson function on a log x scale
    poisson.c <- function (x) {
      (10^x)^(0.5)/(10^x) + min(exprs_cv[gene_valid])
    }

    return(
      ggplot(data.frame(means = log10(expr_mean[gene_valid]),
                        cvs = exprs_cv[gene_valid],
                        gene_type = rep(1, length(expr_mean[gene_valid])) ),
             aes(x = means, y = cvs, col = as.factor(gene_type)) )  +
        geom_point(size = 2, alpha = 0.5) +
        stat_function(fun = poisson.c, col= "red")  +
        scale_colour_manual(values = cbPalette) +
        labs(x = "log10 average molecule count",
             y = ylab,
             title = main)
    )
  }




#' Coefficient of variation versus mean including ERCC genes
plot_poisson_cv <- function(molecules_ENSG, molecules_ERCC,

                            is_log2count = FALSE,
                            include_observed_ERCC = TRUE,
                            main){
  cbPalette <- c("#999999", "#0000FF", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  library(matrixStats)
  if (is_log2count == FALSE) {
    molecules_ENSG <- as.matrix(molecules_ENSG)
    molecules_ERCC <- as.matrix(molecules_ERCC)
    # Remove genes with zero mean molecule count
    which_ENSG_finite <- which(rowMeans(molecules_ENSG) > 0)
    molecules_ENSG <- molecules_ENSG[which_ENSG_finite, ]
    which_ERCC_finite <- which(rowMeans(molecules_ERCC) > 0)
    molecules_ERCC <- molecules_ERCC[which_ERCC_finite, ]

  }
  if (is_log2count == TRUE) {
    molecules_ENSG <- 2^as.matrix(molecules_ENSG)
    molecules_ERCC <- 2^as.matrix(molecules_ERCC)
  }

  # defnine poisson function on a log x scale
  ensg_cv   <- sqrt(rowVars(molecules_ENSG))/rowMeans(molecules_ENSG)
  poisson.c <- function (x) {
    (10^x)^(0.5)/(10^x) + min(ensg_cv)
  }

  # compute the lossy factor based on ERCC
  ####   use maximum likelihood estimate
  ####   dont use the points from ERCC.mol.mean < 0.1 to fit.
  ercc_mean <- rowMeans(molecules_ERCC)
  ercc_cv   <- sqrt(rowVars(molecules_ERCC))/rowMeans(molecules_ERCC)

  require(MASS)
  glm_fit <- glm.nb(round(ercc_mean[log10(ercc_mean) > 0]) ~ 1)
  dispersion <- summary.glm(glm_fit)$dispersion

  # ERCC poisson
  lossy.posson <- function (x) {
    1/sqrt((10^x)/dispersion) + min(ercc_cv)
  }

  # 3 s.d.
  large_sd <- function (x) {
    1/sqrt((10^x)/dispersion/3) + min(ensg_cv)
  }

  if (include_observed_ERCC == TRUE) {
    return(
      ggplot(rbind(data.frame(means = log10(rowMeans(molecules_ENSG)),
                              cvs = sqrt(rowVars(molecules_ENSG))/rowMeans(molecules_ENSG),
                              gene_type = rep(1, NROW(molecules_ENSG) ) ),
                   data.frame(means = log10(ercc_mean),
                              cvs = ercc_cv,
                              gene_type = rep(2, NROW(molecules_ERCC) ) ) ),
             aes(x = means, y = cvs,
                 col = as.factor(gene_type) ) )  +
        geom_point(size = 2, alpha = 0.5) +
        stat_function(fun = poisson.c, col= "red")  +
        stat_function(fun = lossy.posson, col= "blue") +
        stat_function(fun = large_sd, col = "yellow") +
        scale_colour_manual(values = cbPalette) +
        labs(x = "log10 average molecule count",
             y ="Coefficient of variation (CV)",
             title = main)
    )
  }
  if (include_observed_ERCC == FALSE) {
    return(
      ggplot(rbind(data.frame(means = log10(rowMeans(molecules_ENSG)),
                              cvs = sqrt(rowVars(molecules_ENSG))/rowMeans(molecules_ENSG),
                              gene_type = rep(1, NROW(molecules_ENSG))),
                   data.frame(means = log10(ercc_mean),
                              cvs = ercc_cv,
                              gene_type = rep(2, NROW(molecules_ERCC)))),
             aes(x = means, y = cvs, col = as.factor(gene_type)) )  +
        geom_point(size = 2, alpha = 0.5) +
        stat_function(fun = poisson.c, col= "red")  +
        scale_colour_manual(values = cbPalette) +
        labs(x = "log10 average molecule count",
             y ="Coefficient of variation (CV)",
             title = main)
    )
  }
}
