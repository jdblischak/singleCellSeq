#!/bin/bash
set -e

# dedup_umi.py requires pysam.

FILE=$1
BASE=`basename ${FILE%.bam}`
OUTDIR=bam-rmdup-umi

mkdir -p $OUTDIR

if [ ! -s $FILE ]
then
  echo "File is missing or empty: $FILE"
  exit 65
fi

if [ -s $OUTDIR/$BASE.rmdup.bam ]
then
  echo "Output file already exists: $OUTDIR/$BASE.rmdup.bam"
  exit 64
fi

echo "Removing reads with duplicate UMIs..."
dedup_umi.py --method="directional-adjacency" --edit-distance-threshold=1 \
  -I $FILE -v 0 -S $OUTDIR/$BASE.rmdup.bam

echo "Counting number of reads..."
samtools view -c $OUTDIR/$BASE.rmdup.bam > $OUTDIR/$BASE.rmdup.count.txt
