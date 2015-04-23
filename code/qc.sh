#!/bin/bash
set -e

FILE=$1
BASE=`basename ${FILE%.fastq.gz}`
OUTDIR=seqqs

mkdir -p $OUTDIR

if [ ! -s $FILE ]
then
  echo "File is missing or empty: $FILE"
  exit
fi

if [ -s $OUTDIR/${BASE}_len.txt ]
then
  echo "Output file already exists: $OUTDIR/${BASE}_len.txt"
  exit 64
fi

# Run seqqs
zcat $FILE | seqqs -p $OUTDIR/$BASE -
