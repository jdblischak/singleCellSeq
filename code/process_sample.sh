#!/bin/bash
set -e

# Process sequencing data.

# Usage: bash process_sample.sh SOURCE DEST REF_GENOME EXONS

# SOURCE: The directory which contains all the fastq.gz files for a
#   given sample
# DEST: The directory to save the processed data.
# REF_GENOME: The indexed genome for mapping (no extension)
# EXONS: The SAF (simplified annotation format) file for counting reads per exon

# To submit in parallel:
# for DIR in /rawdata/Illumina_Runs/150116_SN_0795_0416_AC5V7FACXX/Demultiplexed/Unaligned/Project_N/*
# do
#   echo "bash process_sample.sh $DIR /mnt/gluster/data/internal_supp/singleCellSeq /mnt/gluster/data/internal_supp/singleCellSeq/genome/combined /mnt/gluster/home/jdblischak/singleCellSeq/data/exons.saf"" | qsub -l h_vmem=12g -V -j y -cwd -N `basename $DIR`
# done

SOURCE=$1
DEST=$2
REF_GENOME=$3 
EXONS=$4 

mkdir -p $DEST/fastq $DEST/bam $DEST/counts

BASE=`basename ${SOURCE}`

echo "Combine reads..."
zcat $SOURCE/*fastq.gz > $DEST/fastq/$BASE.fastq

echo "Run FastQC..."
fastqc $DEST/fastq/$BASE.fastq

echo "Trim UMI..."
umitools trim --end 5 $DEST/fastq/$BASE.fastq NNNNNGGG --verbose > $DEST/fastq/$BASE.trim.fastq

echo "Map with subread-align to genome that combines hg19 and ERCC controls..."
subread-align -uH -i $REF_GENOME -r $DEST/fastq/$BASE.trim.fastq --BAMoutput > $DEST/bam/$BASE.bam

echo "Sort and index bam file..."
samtools sort $DEST/bam/$BASE.bam $DEST/bam/$BASE.sorted
samtools index $DEST/bam/$BASE.sorted.bam

echo "Remove reads with duplicate UMI..."
umitools rmdup $DEST/bam/$BASE.sorted.bam $DEST/bam/$BASE.rmdup.bam > $DEST/bam/$BASE.rmdup.bed

echo "Count reads per gene..."
featureCounts -F SAF -a $EXONS -o $DEST/counts/$BASE.counts.txt $DEST/bam/$BASE.sorted.bam
featureCounts -F SAF -a $EXONS -o $DEST/counts/$BASE.rmdup.counts.txt $DEST/bam/$BASE.rmdup.bam

echo "Remove intermediate files to save space..."
rm $DEST/fastq/$BASE.*fastq $DEST/bam/$BASE.bam
