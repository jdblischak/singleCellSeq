#!/bin/bash

# Runs each single cell through `kallisto pseudoalign`.
# Run from data directory.
# Currently can only run on the same compute node on which it was compiled, spudling07.

# From the README of the modified kallisto
# https://github.com/govinda-kamath/clustering_on_transcript_compatibility_counts/blob/master/modified-kallisto-source/README.md

# Note that we require (path/to/output_file) and not
# (path/to/output_folder). The output file contains the
# transcript-compatibility counts corresponding to the reads that were
# provided in the input. In each line of the output kallisto pseudoalign
# reports the equivalence class id followed by the corresponding
# counts. In the case a read generates a new equivalence class that is
# not in the index, the equivalence class id is set to be the list of
# transcript ids in the new class. This is necessary so that one can
# keep track of new eq.classes generated from different cells, and
# generate the transcript-compatibility counts matrix.

mkdir -p ~/log/run-kallisto.sh

for IND in 19098 19101 19239
do
  for REP in 1 2 3
  do
    for ROW in A B C D E F G H
    do
	for COL in 01 02 03 04 05 06 07 08 09 10 11 12
	do
	  echo "Processing sample $IND.$REP.$ROW$COL"
	  CMD="kallisto pseudoalign -i kallisto/combined.idx \
                               -o kallisto/NA$IND.r$REP.$ROW$COL.txt \
                                trim/$IND.$REP.$ROW$COL*fastq.gz"
          echo $CMD | qsub -l h_vmem=4g -cwd -V -N kallisto.$IND.$REP.$ROW$COL -j y -o ~/log/run-kallisto.sh -l 'hostname=spudling07'
	done
    done
  done
done
