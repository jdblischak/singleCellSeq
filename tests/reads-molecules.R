
library(testit)

assert(
  "file with read counts exists",
  file.exists("../data/reads.txt")
)

assert(
  "file with molecule counts exists",
  file.exists("../data/molecules.txt")
)

reads <- read.table("../data/reads.txt", header = TRUE,
                    stringsAsFactors = FALSE)
molecules <- read.table("../data/molecules.txt", header = TRUE,
                        stringsAsFactors = FALSE)

assert(
  "reads and molecules data have same dimensions",
  dim(reads) == dim(molecules)
)

assert(
  "reads and molecules data have the same column names",
  colnames(reads) == colnames(molecules)
)

assert(
  "reads and molecules data have the same row names",
  rownames(reads) == rownames(molecules)
)

assert(
  "Number of reads per sample greater than or equal to the number of molecules per sample",
  colSums(reads) >= colSums(molecules)
)
