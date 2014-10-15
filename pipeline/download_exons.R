# Download exons from Ensembl for use with Subread featureCounts.

# Currently this includes all genes, which includes non-coding genes.

library("biomaRt")
ensembl <- useMart("ensembl", data = "hsapiens_gene_ensembl")
# attributePages(ensembl)
# atts <- listAttributes(ensembl, page = "feature_page")
# atts[grep("rand", atts$description), ]
exons_all <- getBM(attributes = c("ensembl_gene_id", "ensembl_exon_id",
                                  "chromosome_name", "start_position",
                                  "end_position", "strand"),
                   mart = ensembl)
exons_final <- exons_all[exons_all$chromosome_name %in% c(1:22, "X", "Y", "MT"),
                         c("ensembl_gene_id", "chromosome_name", "start_position",
                           "end_position", "strand")]
colnames(exons_final) <- c("GeneID", "Chr", "Start", "End", "Strand")
# Fix chromosome names
exons_final$Chr <- paste0("chr", exons_final$Chr)
exons_final$Chr <- sub("chrMT", "chrM", exons_final$Chr)
# Fix strand
exons_final$Strand <- ifelse(exons_final$Strand == 1, "+", "-")
# Save as tab-separated file in Simplified Annotation Format (SAF)
write.table(exons_final, "", quote = FALSE, sep = "\t",
            row.names = FALSE)
