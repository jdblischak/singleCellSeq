
library(testit)
suppressPackageStartupMessages(library(dplyr))

assert(
  "file with read counts exists",
  file.exists("../data/reads.txt"),
  file.exists("../data/reads-lcl.txt")
)

assert(
  "file with molecule counts exists",
  file.exists("../data/molecules.txt"),
  file.exists("../data/molecules-lcl.txt")
)

reads <- read.table("../data/reads.txt", header = TRUE,
                    stringsAsFactors = FALSE)
molecules <- read.table("../data/molecules.txt", header = TRUE,
                        stringsAsFactors = FALSE)

reads_lcl <- read.table("../data/reads-lcl.txt", header = TRUE,
                    stringsAsFactors = FALSE)
molecules_lcl <- read.table("../data/molecules-lcl.txt", header = TRUE,
                        stringsAsFactors = FALSE)

assert(
  "reads and molecules data have same dimensions",
  dim(reads) == dim(molecules),
  dim(reads_lcl) == dim(molecules_lcl)
)

assert(
  "reads and molecules data have the same column names",
  colnames(reads) == colnames(molecules),
  colnames(reads_lcl) == colnames(molecules_lcl)
)

assert(
  "reads and molecules data have the same row names",
  rownames(reads) == rownames(molecules),
  rownames(reads_lcl) == rownames(molecules_lcl)
)

assert(
  "iPSCs and LCLs have the same genes measured in the same order",
  rownames(reads) == rownames(reads_lcl),
  rownames(molecules) == rownames(molecules_lcl)
)

assert(
  "Number of reads per sample greater than or equal to the number of molecules per sample",
  colSums(reads) >= colSums(molecules),
  colSums(reads_lcl) >= colSums(molecules_lcl)
)

summary_counts <- read.table("../data/summary-counts.txt", header = TRUE,
                             stringsAsFactors = FALSE)
summary_counts <- arrange(summary_counts, individual, batch, well)

summary_counts_lcl <- read.table("../data/summary-counts-lcl.txt", header = TRUE,
                                 stringsAsFactors = FALSE)
summary_counts_lcl <- arrange(summary_counts_lcl, well)

assert(
  "column sums match Assigned from summary counts",
  colSums(reads) == summary_counts[summary_counts$rmdup == "reads", "Assigned"],
  colSums(molecules) == summary_counts[summary_counts$rmdup == "molecules", "Assigned"],
  colSums(reads_lcl) == summary_counts_lcl[summary_counts_lcl$rmdup == "reads", "Assigned"],
  colSums(molecules_lcl) == summary_counts_lcl[summary_counts_lcl$rmdup == "molecules", "Assigned"]
)
