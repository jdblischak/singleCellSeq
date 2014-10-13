#!/bin/bash

# Create symlinks to fastq files in /rawdata

mkdir -p data/seqs

LANE1=/rawdata/Illumina_Runs/140923_SN_0795_0387_AHA9RUADXX/Unaligned4/Undetermined_indices/Sample_lane1
LANE2=/rawdata/Illumina_Runs/140923_SN_0795_0387_AHA9RUADXX/Unaligned4/Undetermined_indices/Sample_lane2

for fq in $LANE1/*fastq.gz $LANE2/*fastq.gz
do
  #echo $fq
  ln -s $fq /mnt/gluster/home/jdblischak/single-cell-seq/data/seqs/`basename $fq`
done
