#!/bin/bash

#SBATCH -N 1
#SBATCH --mem=5g
#SBATCH -n 1
#SBATCH -t 6-

module load r/4.1.0
Rscript 'C01-ML-Boruta.R' $1
