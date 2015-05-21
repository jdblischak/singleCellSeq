#!/bin/bash
set -e

# Run md5sum on the input file and save in subdirectory "md5sum".

FILE=$1
BASE=`basename ${FILE%.fastq.gz}`
OUTDIR=md5sum

mkdir -p $OUTDIR

if [ ! -s $FILE ]
then
  echo "File is missing or empty: $FILE"
  exit
fi

if [ -s $OUTDIR/$BASE.md5.txt ]
then
  echo "Output file already exists: $OUTDIR/$BASE.md5.txt"
  exit 64
fi

md5sum $FILE > $OUTDIR/$BASE.md5.txt
