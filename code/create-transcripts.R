# Create files with genomic features for examining sequencing coverage.
#
# Usage:
#   Rscript create-transripts.R
#
# Output:
# + tss.bed - Transcription start site per gene
#           - 5' most for forward strand or 3' most for reverse strand
# + tes.bed - Transcription end site per gene
#           - 3' most for forward strand or 5' most for reverse strand
# + transcripts.bed - The most extreme 5' and 3' positions for each gene.
#                   - Corresponds to the tss and tes defined above
#
# Notes:
# + This includes only coding genes, i.e. gene_biotype == "protein_coding"
# + Uses grch37 Ensembl archive
# + Output is in BED format

# Download human transcripts (grch37, hg19)
library("biomaRt")
ensembl <- useMart(host = "grch37.ensembl.org",
                   biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "hsapiens_gene_ensembl")
# attributePages(ensembl)
# [1] "feature_page" "structure" "homologs" "sequences" "snp" "snp_somatic"
atts <- listAttributes(ensembl, page = "feature_page")
# atts[grep("strand", atts$description, ignore.case = TRUE), ]
atts_struct <- listAttributes(ensembl, page = "structure")
# atts_struct[grep("transcript", atts_struct$description, ignore.case = TRUE), ]
transcripts_all <- getBM(attributes = c("ensembl_gene_id", "chromosome_name",
                                  "start_position", "end_position",
                                  "transcript_count", "strand",
                                  "external_gene_name", "gene_biotype"),
                   mart = ensembl)
transcripts <- transcripts_all[transcripts_all$chromosome_name %in%
                                 c(1:22, "X", "Y", "MT") &
                               transcripts_all$gene_biotype == "protein_coding", ]
stopifnot(transcripts$strand %in% c(-1, 1))

# transcript_start and transcript_end are the transcription start and
# transcription end sites for genes on the forward strand, respectively (it is
# the opposite for genes on the reverse strand). start_position is the most 5'
# transcript_start of all the transcripts of a gene. end_position is the 3' most
# transcript_end of all the transcripts of a gene. Thus transcript_start and
# transcript_end can be different for each transcript of a gene, but they always
# have the same start_position and end_position.
# http://stackoverflow.com/a/25434411/2483477

tss <- transcripts
tes <- transcripts
window_size <- 1000


for (g in 1:nrow(transcripts)) {
  if (transcripts$strand[g] == 1) {
    tss_location <- transcripts$start_position[g]
    tes_location <- transcripts$end_position[g]
  } else {
    tss_location <- transcripts$end_position[g]
    tes_location <- transcripts$start_position[g]
  }
  tss$start_position[g] <- tss_location - window_size
  tss$end_position[g] <- tss_location + window_size
  tes$start_position[g] <- tes_location - window_size
  tes$end_position[g] <- tes_location + window_size
}


# Convert to bed file
tss_bed <- data.frame(chr = paste0("chr", tss$chromosome_name),
                      start = tss$start_position - 1,
                      end = tss$end_position,
                      name = tss$ensembl_gene_id,
                      score = tss$transcript_count,
                      strand = ifelse(tss$strand == 1, "+", "-"))
tes_bed <- data.frame(chr = paste0("chr", tes$chromosome_name),
                      start = tes$start_position - 1,
                      end = tes$end_position,
                      name = tes$ensembl_gene_id,
                      score = tes$transcript_count,
                      strand = ifelse(tes$strand == 1, "+", "-"))

transcripts_bed <- data.frame(chr = paste0("chr", transcripts$chromosome_name),
                      start = transcripts$start_position - 1,
                      end = transcripts$end_position,
                      name = transcripts$ensembl_gene_id,
                      score = transcripts$transcript_count,
                      strand = ifelse(transcripts$strand == 1, "+", "-"))


# Fix mitochondrial chromosome name
tss_bed$chr <- sub("chrMT", "chrM", tss_bed$chr)
tes_bed$chr <- sub("chrMT", "chrM", tes_bed$chr)
transcripts_bed$chr <- sub("chrMT", "chrM", transcripts_bed$chr)

# Output
write.table(tss_bed, "../data/tss.bed", quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)
write.table(tes_bed, "../data/tes.bed", quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)
write.table(transcripts_bed, "../data/transcripts.bed", quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)
