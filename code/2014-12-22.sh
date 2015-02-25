#!/bin/bash

FLOW_CELL="/rawdata/Illumina_Runs/141222_SN_0795_0412_AHB79JADXX/Unaligned2"
SUBREAD="/mnt/lustre/home/jdblischak/src/subread-1.4.4-Linux-x86_64/bin"
REF_GENOME="data/genome/combined"
SEQ="data/2014-12-22"
FASTQC="/mnt/lustre/data/tools/FastQC"

# Combine unaligned reads
zcat $FLOW_CELL/Undetermined_indices/Sample_lane[1-2]/*fastq.gz > $SEQ.fastq
gzip $SEQ.fastq

# Map with subread-align to genome that combines hg19, phiX, and ERCC controls
$SUBREAD/subread-align -i $REF_GENOME -r $SEQ.fastq.gz --gzFASTQinput --BAMoutput > $SEQ.bam

# Count reads per feature
samtools view $SEQ.bam | cut -f3 | sort | uniq -c > $SEQ.counts.txt

# Run FastQC
$FASTQC/fastqc -f bam $SEQ.bam
