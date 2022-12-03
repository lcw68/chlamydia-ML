#!/bin/sh

for j in {1..100}; do 
      sbatch -t 7- C01-ML-Boruta-run.sh $j;
done; 