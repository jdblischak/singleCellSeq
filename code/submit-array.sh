#!/bin/bash
set -e

ANALYSIS=$1
MEM=$2
FILES=${*:3}

NAME=`basename $ANALYSIS`
mkdir -p ~/log/$NAME

let numfiles=$#-2
echo $numfiles

echo "#! /bin/bash

files=( ${FILES[@]} )
let index=\$SGE_TASK_ID-1
F=\${files[\$index]}

$ANALYSIS \$F

EXIT_CODE=\$?

if [ \$EXIT_CODE == 0 ]
then
  echo -e \"success\t\$F\"
elif [ \$EXIT_CODE == 64 ]
then
  echo -e \"Output file already exists. Input file was:\t\$F\"
else 
  echo -e \"failure\t\$F\"
fi
" | qsub -l h_vmem=$MEM -V -j y -N $NAME -o ~/log/$NAME -cwd -t 1-$numfiles
