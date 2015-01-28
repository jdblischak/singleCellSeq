#!/bin/bash

# Process sequencing data from run 150116_SN_0795_0416_AC5V7FACXX.

# Usage: bash 2015-01-16.sh SOURCE DEST

# SOURCE: The directory which contains all the fastq.gz files for a
#   given sample
# DEST: The directory to save the processed data.

# To submit in parallel:
# for DIR in /rawdata/Illumina_Runs/150116_SN_0795_0416_AC5V7FACXX/Demultiplexed/Unaligned/Project_N/*
# do
#   echo "bash 2015-01-16.sh $DIR /mnt/gluster/home/ptung/2015-01-16/data" | qsub -l h_vmem=12g -V -j y -cwd -N `basename $DIR`
# done

SOURCE=$1
DEST=$2

# Software
FASTQC="/mnt/lustre/data/tools/FastQC"
UMITOOLS="/mnt/lustre/home/jdblischak/programs/virtualenv-1.11/py2.7/bin/"
SUBREAD="/mnt/lustre/home/jdblischak/src/subread-1.4.4-Linux-x86_64/bin"
SAMTOOLS="/usr/local/bin/"

# Data
REF_GENOME="/mnt/gluster/home/jdblischak/singleCellSeq/pipeline/data/genome/combined"
EXONS="/mnt/gluster/home/jdblischak/singleCellSeq/pipeline/data/genome/exons.saf"

NEW_NAME=$DEST/`basename ${SOURCE}`

echo "Combine reads..."
zcat $SOURCE/*fastq.gz > $NEW_NAME.fastq

echo "Run FastQC..."
$FASTQC/fastqc $NEW_NAME.fastq

echo "Trim UMI..."
$UMITOOLS/umitools trim --end 5 $NEW_NAME.fastq NNNNNGGG --verbose > $NEW_NAME.trim.fastq

echo "Map with subread-align to genome that combines hg19, phiX, and ERCC controls..."
$SUBREAD/subread-align -i $REF_GENOME -r $NEW_NAME.trim.fastq --BAMoutput > $NEW_NAME.bam

echo "Sort and index bam file..."
$SAMTOOLS/samtools sort $NEW_NAME.bam $NEW_NAME.sorted
$SAMTOOLS/samtools index $NEW_NAME.sorted.bam

echo "Remove reads with duplicate UMI..."
$UMITOOLS/umitools rmdup $NEW_NAME.sorted.bam $NEW_NAME.rmdup.bam > $NEW_NAME.rmdup.bed

echo "Count reads per gene..."
$SUBREAD/featureCounts -F SAF -a $EXONS -o $NEW_NAME.counts.txt $NEW_NAME.sorted.bam
$SUBREAD/featureCounts -F SAF -a $EXONS -o $NEW_NAME.rmdup.counts.txt $NEW_NAME.rmdup.bam

echo "Remove intermediate files to save space..."
rm $NEW_NAME.fastq $NEW_NAME.trim.fastq $NEW_NAME.bam
