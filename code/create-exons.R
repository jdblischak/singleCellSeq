# Create exons file for for mapping reads to genes with Subread featureCounts.

# Usage:
#   Rscript create-exons.R ERCC > out.saf
#   Ex:
#   Rscript create-exons.R ../data/ERCC92.gtf > ../data/exons.saf
#
#  ERCC: gtf file of ERCC spike-in controls downloaded from Invitrogen

# Notes:
# + This includes only coding genes, i.e. gene_biotype == "protein_coding"
# + Uses grch37 Ensembl archive
# + Output is in Simplified Annotation Format (SAF)
#     + Columns: GeneID, Chr, Start, End, Strand
#     + Coordinates are 1-based, inclusive on both ends
# + Contains duplicate and overlapping exons (featureCounts handles this)
# + Includes ERRC spike-in controls downloaded from Invitrogen's site
#     + URL: http://media.invitrogen.com.edgesuite.net/softwares/ERCC92.gtf

# Get path to ERCC file
ercc <- commandArgs(trailingOnly = TRUE)

# Download human exons (grch37, hg19)
library("biomaRt")
ensembl <- useMart(host = "grch37.ensembl.org",
                   biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "hsapiens_gene_ensembl")
# attributePages(ensembl)
# [1] "feature_page" "structure" "homologs" "sequences" "snp" "snp_somatic"
atts <- listAttributes(ensembl, page = "feature_page")
# atts[grep("strand", atts$description, ignore.case = TRUE), ]
atts_struct <- listAttributes(ensembl, page = "structure")
# atts_struct[grep("exon", atts_struct$description, ignore.case = TRUE), ]
exons_all <- getBM(attributes = c("ensembl_gene_id", "ensembl_exon_id",
                                  "chromosome_name", "exon_chrom_start",
                                  "exon_chrom_end", "strand",
                                  "external_gene_name",
                                  "gene_biotype"),
                   mart = ensembl)
exons_final <- exons_all[exons_all$chromosome_name %in% c(1:22, "X", "Y", "MT") &
                         exons_all$gene_biotype == "protein_coding",
                         c("ensembl_gene_id", "chromosome_name", "exon_chrom_start",
                           "exon_chrom_end", "strand", "external_gene_name")]
colnames(exons_final) <- c("GeneID", "Chr", "Start", "End", "Strand", "Name")
# Sort by chromosome and position
exons_final <- exons_final[order(exons_final$Chr,
                                 exons_final$Start,
                                 exons_final$End), ]
# Fix chromosome names
exons_final$Chr <- paste0("chr", exons_final$Chr)
exons_final$Chr <- sub("chrMT", "chrM", exons_final$Chr)
# Fix strand
exons_final$Strand <- ifelse(exons_final$Strand == 1, "+", "-")

# Import ERCC data
ercc_gtf <- read.table(ercc, sep = "\t", stringsAsFactors = FALSE)
# http://www.genome.ucsc.edu/FAQ/FAQformat.html#format3
colnames(ercc_gtf) <- c("seqname", # The name of the sequence. Must be a chromosome or scaffold.
                        "source",  # The program that generated this feature.
                        "feature", # The name of this type of feature. Some examples of standard feature types are "CDS", "start_codon", "stop_codon", and "exon".
                        "start",   # The starting position of the feature in the sequence. The first base is numbered 1.
                        "end",     # The ending position of the feature (inclusive).
                        "score",   # A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). If there is no score value, enter ".".
                        "strand",  # Valid entries include '+', '-', or '.' (for don't know/don't care).
                        "frame",   # If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be '.'.
                        "group"    # All lines with the same group are linked together into a single item.
)
ercc_saf <- ercc_gtf[, c("seqname", "seqname", "start", "end", "strand",
                         "seqname")]
colnames(ercc_saf) <- c("GeneID", "Chr", "Start", "End", "Strand", "Name")

# Save as tab-separated file in Simplified Annotation Format (SAF)
saf <- rbind(exons_final, ercc_saf)
write.table(saf, "", quote = FALSE, sep = "\t",
            row.names = FALSE)
