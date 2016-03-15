anno_filter <- read.table("../data/annotation-filter.txt",
                          header = TRUE,
                          stringsAsFactors = FALSE)

molecules_filter <- read.table("../data/molecules-filter.txt",
                               header = TRUE, stringsAsFactors = FALSE)

molecules_final <- read.table("../data/molecules-final.txt",
                              header = TRUE, stringsAsFactors = FALSE)

load("rda/cv-adjusted-summary-pois-expressed/cv.rda")
load("../data/expressed_cv.rda")

gene_ind <- which(rownames(ENSG_cv[[1]]) %in% rownames(expressed_cv[[1]]))
molecules_filter_df <- molecules_filter[which(rownames(molecules_filter) %in%
                                                rownames(expressed_cv[[1]])), ]
molecules_filter_df[which(molecules_filter_df == 0, arr.ind= TRUE)] <- NA

molecules_final_df <- molecules_final[gene_ind, ]
molecules_expressed_df <- molecules_filter_df
molecules_expressed_df[which(molecules_filter_df>0, arr.ind = TRUE)] <- 1


library(Humanzee)
# permuted_data <- Humanzee::permute_cells(log2counts = molecules_final_df,
#                                          grouping_vector = anno_filter$individual,
#                                          number_permute = 4,
#                                          subset_matrix = molecules_expressed_df)
#
# perm_cv <- Humanzee::compute_cv(log2counts = permuted_data[[4]],
#                                 grouping_vector = anno_filter$individual)
#
# perm_cv_adj <- normalize_cv(group_cv = perm_cv,
#                             log2counts = permuted_data[[2]],
#                             anno = anno_filter)
# lapply(perm_cv_adj, function(x) summary(x[[5]]) )
#

#debug(Humanzee::permute_cv_test)
res <- Humanzee::permute_cv_test(log2counts = molecules_final_df,
                subset_matrix = molecules_expressed_df,
                grouping_vector = anno_filter$individual,
                anno = anno_filter,
                number_permute = 4,
                output_rda = FALSE,
                do_parallel = FALSE)
#undebug(Humanzee::permute_cv_test)



# a <- rollapply( data_cv^2[order_gene],
#            width = 30, by = 15,
#            FUN = median, fill = list("extend", "extend", "NA") )
