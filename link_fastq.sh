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

# The sequencing data for the population of cells.
# $ grep 19239 /data/share/PrimaryData/HapGE/RNASeq/list_lanes_formapping_freeze4 
# argonne	090312_HGAC_S100054	2	NA19239	3.5	2
# yale	090121_YOAV_FC310E0_PIPELINE_V2	3	NA19239	2.5	2
# Unfotunately the argonne data does not have the raw sequences, just the
# mapped reads (in some unknown binary format).
ln -s /mnt/lustre/data/share/PrimaryData/HapGE/RNASeq/yale/090121_YOAV_FC310E0_PIPELINE_V2/s_3_sequence.txt.gz data/population/19239_yale.fastq.gz
