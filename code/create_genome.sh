#!/bin/bash
set -e

# Create an indexed genome for mapping sequences derived from human
# mRNA and ERCC spike-in controls.

# + Download sequences (fasta format) for Human (hg19) and ERCC spike-in
# + Build index for mapping with Subread

# First argument is output directory
OUT_DIR=$1

mkdir -p $OUT_DIR
cd $OUT_DIR

# echo "Download hg19 autosomes"
# wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr{1..22}.fa.gz
# 
# echo "Download hg19 X, Y, and M"
# wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr{X,Y,M}.fa.gz
# 
# echo "Unzip fasta files"
# gunzip chr*fa.gz
# 
# echo "Download ERCC data"
# wget http://tools.invitrogen.com/downloads/ERCC92.fa

echo "Build index"
/mnt/gluster/data/internal_supp/singleCellSeq/genome/STAR/bin/Linux_x86_64_static/STAR \
  --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir star_index \
  --sjdbGTFfile hg19-ensembl-genes.gtf \
  --genomeFastaFiles *fa \
  --sjdbOverhang 99
