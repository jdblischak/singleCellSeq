#!/bin/bash

FASTQGZ=$1

SUBREAD=/mnt/lustre/home/jdblischak/src/subread-1.4.4-Linux-x86_64/bin
FASTQC=/mnt/lustre/data/tools/FastQC

echo "Unzip fastq file"
FASTQ=${FASTQGZ%.gz}
zcat $FASTQGZ > $FASTQ

echo "Assess quality with FastQC"
$FASTQC/fastqc $FASTQ

echo "Trim UMI's from 5' end of reads"
umitools trim --end 5 $FASTQ NNNNNGGG --verbose > $FASTQ.trim

echo "Map reads"
BAM=${FASTQ%fastq}bam
$SUBREAD/subread-align -i data/genome/combined -r $FASTQ.trim --BAMoutput > $BAM

echo "Assess mapping"
COUNTS=${FASTQ%fastq}counts.txt
samtools view $BAM | cut -f 3 | sort | uniq -c > $COUNTS

echo "Sort and index mapped reads"
samtools sort $BAM ${BAM%bam}sorted
samtools index ${BAM%bam}sorted.bam

echo "Save one read per UMI"
umitools rmdup ${BAM%bam}sorted.bam $BAM.umi > ${BAM%bam}umi.diffs.bed

echo "Remove intermediate file"
rm $FASTQ $FASTQ.trim
