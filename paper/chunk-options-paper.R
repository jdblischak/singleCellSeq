# For more information on available chunk options, see
# http://yihui.name/knitr/options#chunk_options

library("knitr")
opts_chunk$set(include = FALSE)

# Numbering for main figures
fig_main_qc <- "1"
fig_main_subsample <- "2"
fig_main_batch <- "3"
fig_main_normalization <- "4"
fig_main_noise <- "5"

# Numbering for supplementary figures
fig_supp_qc <- "S1"
fig_supp_lda <- "S2"
fig_supp_subsample <- "S3"
fig_supp_variance <- "S4"
fig_supp_dropout <- "S5"
fig_supp_permutation <- "S6"
fig_supp_noisygene <- "S7"
fig_supp_plurigene <- "S8"
fig_supp_design <- "S9"
fig_supp_proportion <- "S10"
fig_supp_CV <- "S11"

# Numbering for supplementary tables
table_collection <- "S1"
table_qualitycell <- "S2"
table_noisygene <- "S3"
table_GO <- "S4"

# Import data
anno_filter <- read.table("../data/annotation-filter.txt", header = TRUE,
                          stringsAsFactors = FALSE)
molecules_filter <- read.table("../data/molecules-filter.txt", header = TRUE,
                               stringsAsFactors = FALSE)
