#!/bin/bash
set -e

FILE=$1
BASE=`basename ${FILE%.bam}`
OUTDIR=bam-processed

mkdir -p $OUTDIR

if [ ! -s $FILE ]
then
  echo "File is missing or empty: $FILE"
  exit 65
fi

if [ -s $OUTDIR/$BASE.sorted.bam ]
then
  echo "Output file already exists: $OUTDIR/$BASE.sorted.bam"
  exit 64
fi

echo "Sorting file..."
samtools sort $FILE $OUTDIR/$BASE.sorted

echo "Indexing file..."
samtools index $OUTDIR/$BASE.sorted.bam

echo "Counting number of reads..."
# Only count mapped reads. Unmapped reads have score of 0
samtools view -c -q 1 $OUTDIR/$BASE.sorted.bam > $OUTDIR/$BASE.processed.count.txt
