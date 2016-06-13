#!/bin/bash
set -e

FILE=$1
BASE=`basename ${FILE%.fastq.gz}`
OUTDIR=sickle

mkdir -p $OUTDIR

if [ ! -s $FILE ]
then
  echo "File is missing or empty: $FILE"
  exit 65
fi

if [ -s $OUTDIR/$BASE.sickle.fastq.gz ]
then
  echo "Output file already exists: $OUTDIR/$BASE.sickle.fastq.gz"
  exit 64
fi

# Run sickle only cutting from the 3' end
sickle se -f <(zcat $FILE) -t sanger -o $OUTDIR/$BASE.sickle.fastq.gz -x -g

bioawk -c fastx 'END{print NR}' $OUTDIR/$BASE.sickle.fastq.gz > $OUTDIR/$BASE.sickle.count.txt
