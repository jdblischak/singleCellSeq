#!/bin/bash
set -e

FILE=$1
BASE=`basename ${FILE%.bam}`
OUTDIR=bam-processed

mkdir -p $OUTDIR

if [ ! -s $FILE ]
then
  echo "File is missing or empty: $FILE"
  exit
fi

if [ -s $OUTDIR/$BASE.sorted.bam ]
then
  echo "Output file already exists: $OUTDIR/$BASE.sorted.bam"
  exit 64
fi

echo "Sorting file..."
samtools sort <(samtools view -b -q 10 $FILE) $OUTDIR/$BASE.sorted

echo "Indexing file..."
samtools index $OUTDIR/$BASE.sorted.bam

echo "Counting number of reads..."
samtools view -c $OUTDIR/$BASE.sorted.bam > $OUTDIR/$BASE.processed.count.txt
