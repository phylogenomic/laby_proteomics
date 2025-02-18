#!/bin/bash
#SBATCH --nodes=1
#SBATCH -p extended-96core
#SBATCH --time=7-00:00:00

# sbatch --export=class='laby' --job-name=proteomes  bin/1.2.Proteome_align.sl
cd /gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/phylogenomics

# Orthofinder
module load phylo/1.0
cd data
orthofinder -f proteomes -M msa -t 96
#cd proteomes/OrthoFinder/Results_MSA
