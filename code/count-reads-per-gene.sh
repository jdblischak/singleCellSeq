#!/bin/bash
set -e

FILE=$1
BASE=`basename ${FILE%.bam}`
EXONS=genome/exons.saf
OUTDIR=counts

mkdir -p $OUTDIR

if [ ! -s $FILE ]
then
  echo "File is missing or empty: $FILE"
  exit 65
fi

if [ -s $OUTDIR/$BASE.genecounts.txt ]
then
  echo "Output file already exists: $OUTDIR/$BASE.genecounts.txt"
  exit 64
fi

echo "Counting reads per gene..."
featureCounts -a $EXONS -F SAF -R -o $OUTDIR/$BASE.genecounts.txt $FILE
