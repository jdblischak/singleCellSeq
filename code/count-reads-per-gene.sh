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
featureCounts -a $EXONS -F SAF -R -s 1 -o $OUTDIR/$BASE.genecounts.txt $FILE

# Bug introduced in 1.5 series. In 1.4 series, the "-R" file was the input
# file with .featureCounts appended. Now it writes to the working directory
# and concatenates the path to the beginning of the filename.
# https://groups.google.com/forum/#!topic/subread/BjVCMxECE3o
# I move it to the output directory since that makes more sense anyways.
mv *$BASE.bam.featureCounts $OUTDIR/$BASE.bam.featureCounts
