library(Humanzee)
library(matrixStats)
load("/home/joycehsiao/scratch-midway/single-cell/permute-cv-express/rda-for-midway-expressed-cells.rda")
load("/home/joycehsiao/scratch-midway/single-cell/permute-cv-express/permuted-distance.rda")

permuted_pval <- data.frame(mad_pval = rowMeans(permuted_distance > mad_expressed) )
save(permuted_pval, 
	file = "/home/joycehsiao/scratch-midway/single-cell/permute-cv-express/permuted-pval-expressed.rda")
