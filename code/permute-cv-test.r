library(Humanzee)
#load("rda/cv-adjusted-statistical-test-permute/permute-cv-test.rda")
load("/home/joycehsiao/scratch-midway/single-cell/permute-cv/permute-cv-test.rda")

permute_cv_expressed(log2counts = molecules_final, 
		  grouping_vector = anno_filter$individual, 
		  anno = anno_filter, 
		  number_permute = dim(molecules_final)[1],
                  output_rda = TRUE,
                  do_parallel = TRUE,
                  number_cores = 12)
