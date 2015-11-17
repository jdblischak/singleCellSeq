#!/bin/bash

#SBATCH --job-name=permute-cv-compute-pval
#SBATCH --output=permute-cv-compute-pval.out
#SBATCH --error=permute-cv-compute-pval.err
#SBATCH --time=36:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joycehsiao@uchicago.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --partition=westmere

module add R/3.2
module load R/3.2

Rscript permute-cv-compute-pval.r

