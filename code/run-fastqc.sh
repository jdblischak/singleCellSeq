#!/bin/bash

FLOW_CELL=$1
SEQ=$2
FASTQC="/mnt/lustre/data/tools/FastQC"

# Combine all reads
FILES=`find $FLOW_CELL -name "*fastq.gz"`
zcat $FILES > $SEQ.fastq

# Run FastQC
$FASTQC/fastqc $SEQ.fastq

rm $SEQ.fastq
