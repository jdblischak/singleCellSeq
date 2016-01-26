library(Humanzee)
library(matrixStats)
#load("/home/joycehsiao/scratch-midway/single-cell/permute-cv/permute-cv-test.rda")
load("/home/joycehsiao/scratch-midway/single-cell/permute-cv/adj-cv.rda")
load("/home/joycehsiao/scratch-midway/single-cell/permute-cv/permuted-distance.rda")

mad_permuted <- as.matrix( do.call(cbind, lapply(permuted_distance, "[[", 1)) )

permuted_pval <- data.frame(mad_pval = rowMeans(mad_permuted > mad))
save(permuted_pval, 
	file = "/home/joycehsiao/scratch-midway/single-cell/permute-cv/permuted-pval.rda")
