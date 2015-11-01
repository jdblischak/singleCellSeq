library(Humanzee)
load("./permute-cv-test.rda")

permute_cv_test(log2counts = molecules_ENSG, 
		  grouping_vector = anno_filter$individual, 
		  anno = anno_filter, 
		  number_permute = 1000,
                  output_rda = TRUE,
                  do_parallel = TRUE,
                  number_cores = 12)
