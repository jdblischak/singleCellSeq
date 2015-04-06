#!/bin/bash

# Run snakemake.

# The first argument decides whether or not to submit the jobs to
# qsub. If it is "t", jobs are submitted. Any other string will result
# in no submission. The remaining arguments are added to the snakemake
# command.

# To submit this script:
# qsub -l h_vmem=8g -cwd -N snakemake -V -j y snakefile.sh t

SUBMIT=$1
ARGS=${*:2}

# Command parts
BASE="snakemake --ri -kps snakefile.py"
CLUSTER="-j 96 -c \"ssh spudhead 'qsub -l h_vmem={params.h_vmem} -l bigio={params.bigio} -N {params.name} -V -j y -wd /mnt/gluster/data/internal_supp/singleCellSeq/ -o {log}'\""
#CLUSTER="-j 96 -c \"qsub -l h_vmem={params.h_vmem} -l bigio={params.bigio} -N {params.name} -V -j y -cwd -o {log}\""
WORKDIR="-d /mnt/gluster/data/internal_supp/singleCellSeq/"

if [ $SUBMIT == "t" ]
then
  CMD="$BASE $WORKDIR $CLUSTER $ARGS"
  echo $CMD
  eval $CMD
else
  CMD="$BASE $WORKDIR $ARGS"
  echo $CMD
  eval $CMD
fi
