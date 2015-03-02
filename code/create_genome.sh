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

echo "Download human genome"
rsync -avzuP globus.opensciencedatacloud.org::public/illumina/igenomes/Homo_sapiens/UCSC/hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa .
mv genome.fa hg19.fa

echo "Download ERCC data"
wget http://tools.invitrogen.com/downloads/ERCC92.fa

echo "Build index"
subread-buildindex -o combined hg19.fa ERCC92.fa

echo "Remove intermediate files"
rm -rf hg19.fa
