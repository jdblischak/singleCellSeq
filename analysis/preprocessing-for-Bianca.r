#' Preprocessing steps for single cell data
#' 
#' This script contains steps we currently take to filter and normalize 
#' single cell data. These steps are our
#' preliminary attemps at processing, and are by no means complete.
#' 
#' Updated: 2015-09-26
#' 
#' Usage: Default Working directory is set to  be "singleCellSeq/analysis."
#'        We recommend cloning the singleCellSeq GitHub repo, and then Working
#'        from your local singleCellSeq directory.


##---- Set up ----#
library(edgeR)


##---- Import data ----#

## Import annotation
anno <- read.table("../data/annotation.txt", header = TRUE,
                   stringsAsFactors = FALSE)

## Import read counts
reads <- read.table("../data/reads.txt", header = TRUE,
                    stringsAsFactors = FALSE)

## Import molecule counts
molecules <- read.table("../data/molecules.txt", header = TRUE, stringsAsFactors = FALSE)


##---- Quality control and filter ----#

## Remove bulk samples
single_samples <- anno$well != "bulk"
anno_single <- anno[ which(single_samples), ]
molecules_single <- molecules[ , which(single_samples)]
reads_single <- reads[ , which(single_samples)]
stopifnot(ncol(molecules_single) == nrow(anno_single),
          colnames(molecules_single) == anno_single$sample_id)

## Remove ERCC genes
ii_nonERCC <- grep("ERCC", rownames(molecules_single), invert = TRUE)
molecules_single_ENSG <- molecules_single[ii_nonERCC, ]


## Remove gene with a sum of 0 count across cells
expressed_single_ENSG <- rowSums(molecules_single_ENSG) > 0
molecules_single_ENSG <- molecules_single_ENSG[expressed_single_ENSG, ]

## Remove gene with molecule count larger than 1024 (15 if them)
overexpressed_genes <- rownames(molecules_single_ENSG)[apply(molecules_single_ENSG, 1,
                                                        function(x) any(x >= 1024))]
molecules_single_ENSG <- molecules_single_ENSG[!(rownames(molecules_single_ENSG) %in% overexpressed_genes), ]

## Collision probability correction (output count data)
molecules_single_cpm <- cpm(molecules_single_collision, log = TRUE)

## CPM normalization of molecule counts
molecules_single_collision <- -1024 * log(1 - molecules_single_ENSG / 1024)


