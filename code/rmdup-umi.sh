#!/bin/bash
set -e

FILE=$1
BASE=`basename ${FILE%.bam}`
OUTDIR=bam-rmdup-umi

mkdir -p $OUTDIR

if [ ! -s $FILE ]
then
  echo "File is missing or empty: $FILE"
  exit
fi

if [ -s $OUTDIR/$BASE.rmdup.bam ]
then
  echo "Output file already exists: $OUTDIR/$BASE.rmdup.bam"
  exit
fi

echo "Removing reads with duplicate UMIs..."
umitools rmdup $FILE $OUTDIR/$BASE.rmdup.bam > $OUTDIR/$BASE.rmdup.bed

echo "Counting number of reads..."
samtools view -c $OUTDIR/$BASE.rmdup.bam > $OUTDIR/$BASE.rmdup.count.txt

echo -e "success\t$BASE"
