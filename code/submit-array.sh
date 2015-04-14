#!/bin/bash
set -e

ANALYSIS=$1
MEM=$2
FILES=${*:3}

mkdir -p log

let numfiles=$#-2
echo $numfiles

echo "#! /bin/bash

files=( ${FILES[@]} )
let index=\$SGE_TASK_ID-1
F=\${files[\$index]}

$ANALYSIS \$F

if [ \$? == 0 ]
then
  echo -e \"success\t\$F\"
elif [ \$? == 64 ]
then
  echo -e \"Output file already exists\t\$F\"
else 
  echo -e \"failure\t\$F\"
fi
" | qsub -l h_vmem=$MEM -V -j y -N `basename $ANALYSIS` -o log -cwd -t 1-$numfiles
