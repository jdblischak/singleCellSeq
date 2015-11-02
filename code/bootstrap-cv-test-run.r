library(Humanzee)
load("./bootstrap-cv-test.rda")

bootstrap_cv_test(log2counts = molecules_ENSG, 
		  grouping_vector = anno_filter$individual, 
		  anno = anno_filter, 
		  num_bootstrap = 1000,
                  output_rda = TRUE,
                  do_parallel = TRUE,
                  number_cores = 12)
