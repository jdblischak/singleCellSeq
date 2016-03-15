#!/bin/bash
set -e

# Create an indexed transcriptome for pseudoaligning sequences derived
# from human mRNA and ERCC spike-in controls.

# + Download transcript sequences (fasta format) for Human (hg38) and ERCC spike-in
# + Build index for pseudoaligning with kallisto

# First argument is output directory
OUT_DIR=$1

mkdir -p $OUT_DIR
cd $OUT_DIR

echo "Download hg38 transcriptome"
wget wget http://bio.math.berkeley.edu/kallisto/transcriptomes/Homo_sapiens.GRCh38.rel79.cdna.all.fa.gz

echo "Download ERCC data"
wget http://tools.invitrogen.com/downloads/ERCC92.fa

echo "Build index"
kallisto index -i combined.idx ERCC92.fa Homo_sapiens.GRCh38.rel79.cdna.all.fa.gz
