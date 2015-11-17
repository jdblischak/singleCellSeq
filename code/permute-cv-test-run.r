library(Humanzee)
load("rda/cv-adjusted-statistical-test-permute/permute-cv-test.rda")

permute_cv_test(log2counts = molecules_ENSG, 
		  grouping_vector = anno_filter$individual, 
		  anno = anno_filter, 
		  number_permute = 10000,
                  output_rda = TRUE,
                  do_parallel = TRUE,
                  number_cores = 12)
