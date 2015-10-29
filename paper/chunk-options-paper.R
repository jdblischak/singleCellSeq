# For more information on available chunk options, see
# http://yihui.name/knitr/options#chunk_options

library("knitr")
opts_chunk$set(include = FALSE)

# Numbering for main figures
fig_main_qc <- "1"
fig_main_cpm <- "2"
fig_main_subsample <- "3"
fig_main_normalization <- "4"
fig_main_noise <- "5"

# Numbering for supplementary figures
fig_supp_flowcell <- "S1"
fig_supp_sickle <- "SX"

# Import data
anno_filter <- read.table("../data/annotation-filter.txt", header = TRUE,
                          stringsAsFactors = FALSE)
molecules_filter <- read.table("../data/molecules-filter.txt", header = TRUE,
                               stringsAsFactors = FALSE)
