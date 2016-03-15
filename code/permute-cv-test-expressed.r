library(Humanzee)
load("/home/joycehsiao/scratch-midway/single-cell/permute-cv-express/rda-for-midway-expressed-cells.rda")

permute_cv_test(log2counts = molecules_final_df, 
			    subset_matrix = molecules_expressed_df,
	  			grouping_vector = anno_filter$individual, 
			    anno = anno_filter, 
#			    number_permute = 10,			    
			    number_permute = dim(molecules_final_df)[1],
                output_rda = TRUE,
                do_parallel = TRUE,
                number_cores = 12)
