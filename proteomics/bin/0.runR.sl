#!/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=long-40core
#SBATCH --time=8:00:00
#SBATCH --ntasks-per-node=40

#RUN WITH:

echo "Start time:" 
date

module load R/4.2.1
#Rscript bin/7.0.CombineAnnotations.R
Rscript bin/7.4.DEP_independentsets.R

echo "End time:" 
date