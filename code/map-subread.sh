#!/bin/bash
set -e

FILE=$1
GENOME=genome/combined
BASE=`basename ${FILE%.fastq.gz}`
OUTDIR=bam

mkdir -p $OUTDIR

if [ ! -s $FILE ]
then
  echo "File is missing or empty: $FILE"
  exit
fi

if [ -s $OUTDIR/$BASE.bam ]
then
  echo "Output file already exists: $OUTDIR/$BASE.bam"
  exit 64
fi

subread-align -i $GENOME -r $FILE --gzFASTQinput --BAMoutput -uH > $OUTDIR/$BASE.bam
