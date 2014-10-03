files <- list.files("data", "counts.txt", full.name = TRUE)
seqs <- scan(files[1], what = "character", quiet = TRUE)[c(FALSE, TRUE)]
seqs <- sub("\\*", "unknown", seqs)
counts <- matrix(nrow = length(seqs), ncol = length(files))
rownames(counts) <- seqs
colnames(counts) <- files
for (i in seq_along(files)) {
  f <- scan(files[i], what = "character", quiet = TRUE)
  counts[, i] <- as.numeric(f[c(TRUE, FALSE)])
}

total_human <- sum(counts[grep("chr", rownames(counts)), ])
total_ercc <- sum(counts[grep("ERCC", rownames(counts)), ])
total_phix <- sum(counts["phix", ])
total_unknown <- sum(counts["unknown", ])

cat(sprintf("Human reads:\t%8d\n", total_human))
cat(sprintf("ERCC reads:\t%8d\n", total_ercc))
cat(sprintf("PhiX reads:\t%8d\n", total_phix))
cat(sprintf("Unknown reads:\t%8d\n", total_unknown))
