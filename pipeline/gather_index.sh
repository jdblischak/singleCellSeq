#!/bin/bash

# Extracts the index sequences from multiple fastq.gz files and writes them
# in fasta format.

FILES=`find /rawdata/Illumina_Runs/141030_SN_0795_0391_AHAV20ADXX/Unaligned -name "*fastq.gz"`
#FILES=/rawdata/Illumina_Runs/141030_SN_0795_0391_AHAV20ADXX/Unaligned/Project_N/Sample_19239_LCL_01/19239_LCL_01_CGTCTAAT_L001_R1_001.fastq.gz
FC=HAV20ADXX
OUT=index.fa

rm $OUT

# For each fastq file, add its index sequences to output file $OUT
for F in $FILES
do
  echo $F
  zcat $F | grep $FC > headers.txt
  cat headers.txt | sed s/@/"> @"/ | sed s/" 1:N:0:"/"\n"/ >> $OUT
done

# Clean up
rm headers.txt
