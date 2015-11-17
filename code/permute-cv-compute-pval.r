library(Humanzee)
library(matrixStats)
load("rda/cv-adjusted-statistical-test-permute/permute-cv-test.rda")
load("rda/cv-adjusted-statistical-test-permute/adj-cv.rda")
load("rda/cv-adjusted-statistical-test-permute/permuted-distance.rda")

squared_dev <- do.call(cbind, lapply(permuted_distance, "[[", 1))
abs_dev <- do.call(cbind, lapply(permuted_distance, "[[", 2))

permuted_pval <- data.frame(squared_dev = rowMeans(squared_dev > df_norm$squared_dev),
			   abs_dev = rowMeans(abs_dev > df_norm$abs_dev) )

save(permuted_pval, 
	file = "rda/cv-adjusted-statistical-test-permute/permuted-pval.rda")
