files <- list.files("data", "counts.txt", full.name = TRUE)
total_human <- 0
total_ercc <- 0
total_phix <- 0
total_unknown <- 0
for (f in files) {
  raw <- scan(f, what = "character", quiet = TRUE)
  counts <- as.numeric(raw[c(TRUE, FALSE)])
  names(counts) <- raw[c(FALSE, TRUE)]
  total_human <- total_human + sum(counts[grep("chr", names(counts))])
  total_ercc <- total_ercc + sum(counts[grep("ERCC", names(counts))])
  total_phix <- total_phix + sum(counts[grep("phix", names(counts))])
  total_unknown <- total_unknown + sum(counts[grep("\\*", names(counts))])
}

stopifnot(total_human + total_ercc + total_phix + total_unknown
          == 4e6 * length(files))

cat(sprintf("Human reads:\t%8d\n", total_human))
cat(sprintf("ERCC reads:\t%8d\n", total_ercc))
cat(sprintf("PhiX reads:\t%8d\n", total_phix))
cat(sprintf("Unknown reads:\t%8d\n", total_unknown))
