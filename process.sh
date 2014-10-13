#!/bin/bash

FASTQGZ=$1

SUBREAD=/mnt/lustre/home/jdblischak/src/subread-1.4.4-Linux-x86_64/bin
FASTQC=/mnt/lustre/data/tools/FastQC

echo "Unzip fastq file"
FASTQ=${FASTQGZ%.gz}
zcat $FASTQGZ > $FASTQ

echo "Assess quality with FastQC"
$FASTQC/fastqc $FASTQ

echo "Map reads"
BAM=${FASTQ%fastq}bam
$SUBREAD/subread-align -i data/genome/combined -r $FASTQ --BAMoutput --trim5 8 > $BAM

echo "Assess mapping"
COUNTS=${FASTQ%fastq}counts.txt
samtools view $BAM | cut -f 3 | sort | uniq -c > $COUNTS

echo "Remove intermediate file"
rm $FASTQ
