library("testit")
library("data.table")
suppressPackageStartupMessages(library("dplyr"))

summary_counts <- fread("../data/summary-counts.txt")
setDF(summary_counts)

# Test basics of featureCounts columns -----------------------------------------

assert(
  "Unused count columns are all zero",
  summary_counts$Unassigned_MultiMapping == 0,
  summary_counts$Unassigned_MappingQuality == 0,
  summary_counts$Unassigned_FragmentLength == 0,
  summary_counts$Unassigned_Chimera == 0,
  summary_counts$Unassigned_Secondary == 0,
  summary_counts$Unassigned_Nonjunction == 0,
  summary_counts$Unassigned_Duplicate == 0
)

assert(
  "Unassigned ambiguous counts always less than mapped counts",
  summary_counts$Unassigned_Ambiguity < summary_counts$Assigned
)

assert(
  "Unassigned unmapped molecule counts are zero because filtered in rmdup UMI step",
  summary_counts$Unassigned_Unmapped[summary_counts$rmdup == "molecules"] == 0
)

# Test concordance with per gene-by-sample counts files

reads <- read.table("../data/reads.txt", header = TRUE)

molecules <- read.table("../data/molecules.txt", header = TRUE)

assert(
  "Sum of read counts are equal to Assigned summary counts",
  ncol(reads) == nrow(summary_counts[summary_counts$rmdup == "reads", ]),
  colSums(reads) == summary_counts$Assigned[summary_counts$rmdup == "reads"]
)

assert(
  "Sum of molceule counts are equal to Assigned summary counts",
  ncol(molecules) == nrow(summary_counts[summary_counts$rmdup == "molecules", ]),
  colSums(molecules) == summary_counts$Assigned[summary_counts$rmdup == "molecules"]
)
