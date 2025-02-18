#!/usr/bin/env bash

# module load slurm
# module load shared
# sbatch download

#SBATCH --job-name=download
#SBATCH --output=diamond_nr.log
#SBATCH --ntasks-per-node=24
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH -p long-24core

module load diamond/2.0.10
module load anaconda/3

source activate myenv
esearch  -db "protein" -query "Stramenopiles[Organism]" | efetch -format fasta > seq_stramenopiles.fa
