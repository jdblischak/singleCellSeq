#!/bin/bash

#SBATCH --job-name=bootstrap-cv-test
#SBATCH --output=bootstrap-cv-test.out
#SBATCH --error=bootstrap-cv-test.err
#SBATCH --time=36:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joycehsiao@uchicago.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --partition=westmere

module add R/3.2
module load R/3.2

Rscript bootstrap-cv-test-run.r

