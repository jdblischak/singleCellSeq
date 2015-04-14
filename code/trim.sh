#!/bin/bash
set -e

FILE=$1
BASE=`basename ${FILE%.fastq.gz}`
OUTDIR=trim

mkdir -p $OUTDIR

if [ ! -s $FILE ]
then
  echo "File is missing or empty: $FILE"
  exit
fi

if [ -s $OUTDIR/$BASE.trim.fastq.gz ]
then
  echo "Output file already exists: $OUTDIR/$BASE.trim.fastq.gz"
  exit 64
fi

umitools trim --end 5 <(zcat $FILE) NNNNNGGG --verbose --top 50 \
  2> $OUTDIR/$BASE.trim.stats.txt | gzip -c 1> $OUTDIR/$BASE.trim.fastq.gz

zcat $OUTDIR/$BASE.trim.fastq.gz | grep "@" | wc -l > $OUTDIR/$BASE.trim.count.txt
