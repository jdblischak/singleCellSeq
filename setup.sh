#!/bin/bash

# + Download sequences (fasta format) for Human (hg19), PhiX, and ERCC spike-in
# + Build index for mapping with Subread

SUBREAD=/mnt/lustre/home/jdblischak/src/subread-1.4.4-Linux-x86_64/bin/

mkdir -p data/genome
cd data/genome

echo "Download human genome"
rsync -avzuP globus.opensciencedatacloud.org::public/illumina/igenomes/Homo_sapiens/UCSC/hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa .
mv genome.fa hg19.fa

echo "Download phiX genome"
rsync -avzuP globus.opensciencedatacloud.org::public/illumina/igenomes/PhiX/NCBI/1993-04-28/PhiX_NCBI_1993-04-28.tar.gz .
tar -xzf PhiX_NCBI_1993-04-28.tar.gz

echo "Download ERCC data"
wget http://tools.invitrogen.com/downloads/ERCC92.fa

echo "Build index"
$SUBREAD/subread-buildindex -o combined hg19.fa PhiX/NCBI/1993-04-28/Sequence/WholeGenomeFasta/genome.fa ERCC92.fa

echo "Remove intermediate files"
rm -r hg19.fa PhiX_NCBI_1993-04-28.tar.gz README.txt PhiX ERCC92.fa
