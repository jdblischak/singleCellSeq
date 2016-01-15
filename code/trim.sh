#!/bin/bash
set -e

FILE=$1
BASE=`basename ${FILE%.fastq.gz}`
OUTDIR=trim

# Directory to send file with reads w/o UMIs
INVALID_DIR=invalid

mkdir -p $OUTDIR $INVALID_DIR

if [ ! -s $FILE ]
then
  echo "File is missing or empty: $FILE"
  exit 65
fi

if [ -s $OUTDIR/$BASE.trim.fastq.gz ]
then
  echo "Output file already exists: $OUTDIR/$BASE.trim.fastq.gz"
  exit 64
fi

umitools trim --end 5 --verbose --top 50 $FILE NNNNNGGG \
  --invalid $INVALID_DIR/$BASE.invalid.fastq \
  2> $INVALID_DIR/$BASE.trim.stats.txt | gzip -c > $OUTDIR/$BASE.trim.fastq.gz

bioawk -c fastx 'END{print NR}' $OUTDIR/$BASE.trim.fastq.gz > $OUTDIR/$BASE.trim.count.txt

bioawk -c fastx 'END{print NR}' $INVALID_DIR/$BASE.invalid.fastq > $INVALID_DIR/$BASE.invalid.count.txt

gzip $INVALID_DIR/$BASE.invalid.fastq
