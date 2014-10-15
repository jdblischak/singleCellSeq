# Combine the feature counts into a count matrix.

# All the counts files end with counts.txt. Read in the first to get the
# annotation.
files <- list.files(c("data/seqs", "data/population"), pattern = "counts.txt$",
                    full.names = TRUE)
gene_anno <- read.table(files[1], header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)
gene_anno <- gene_anno[, c("Geneid", "Chr", "Start", "End", "Strand", "Length")]
gene_counts <- matrix(0, nrow = nrow(gene_anno), ncol = 3)
colnames(gene_counts) <- c("two_cell_reads", "two_cell_umi", "NA19239")

# Read in each file, check that it has the same rows in the correct order,
# then sum if it is part of the two_cell data or simply add to count matrix.
for (fname in files) {
  f <- read.table(fname, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  stopifnot(f$Geneid == gene_anno$Geneid)
  if (grepl("Undetermined", fname)) {
    gene_counts[, "two_cell_reads"] <- gene_counts[, "two_cell_reads"] + f[, 7]
    gene_counts[, "two_cell_umi"] <- gene_counts[, "two_cell_umi"] + f[, 8]
  } else if (grepl("19239", fname)) {
    gene_counts[, "NA19239"] <- f[, 7]
  } else {
    stop("Invalid filename")
  }
}

# Write annotation plus counts to standard output
write.table(cbind(gene_anno, gene_counts), file = "", quote = FALSE, sep = "\t",
            row.names = FALSE)
