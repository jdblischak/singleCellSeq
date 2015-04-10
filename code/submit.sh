#!/bin/bash
set -e

ANALYSIS=$1
MEM=$2
FILES=${*:3}

mkdir -p log

for F in $FILES
do
  qsub -l h_vmem=$MEM -V -j y -N `basename $ANALYSIS`.`basename $F` -o log -cwd -b y $ANALYSIS $F
done

echo success.`basename $ANALYSIS`
