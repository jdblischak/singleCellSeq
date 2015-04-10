#!/bin/bash
set -e

FILE=$1
GENOME=genome/star_index
BASE=`basename ${FILE%.fastq.gz}`
OUTDIR=bam

mkdir -p $OUTDIR

if [ ! -s $FILE ]
then
  echo "File is missing or empty: $FILE"
  exit
fi

if [ -s $OUTDIR/$BASE.Aligned.sortedByCoord.out.bam ]
then
  echo "Output file already exists: $OUTDIR/$BASE.Aligned.sortedByCoord.out.bam"
  exit
fi

STAR \
  --runThreadN 8 \
  --genomeDir $GENOME \
  --readFilesIn	$FILE \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix $OUTDIR/$BASE.

# STAR \
#   --runThreadN 8 \
#   --genomeDir $GENOME \
#   --readFilesIn	$FILE \
#   --readFilesCommand zcat \
#   --outSAMtype BAM SortedByCoordinate \
#   --outFilterScoreMin 255 \
#   --outFileNamePrefix $OUTDIR/$BASE.

#  --genomeLoad LoadAndKeep \
#  --limitBAMsortRAM 10000000000 \

samtools view -c -q 255 $OUTDIR/$BASE.Aligned.sortedByCoord.out.bam > $OUTDIR/$BASE.map.count.txt

echo success.$BASE
