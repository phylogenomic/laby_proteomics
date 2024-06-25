#!/bin/sh

#SBATCH --job-name=multiloc
#SBATCH -o multiloc.sh.out
#SBATCH --nodes=1
#SBATCH -p extended-40core
#SBATCH --time=7-00:00:00

# AUTHOR: ALEX GIL GOMEZ
# INPUTS: GeneCatalog
# OUTPUT: Multiloc file
# INFO: Runs Multiloc on gene catalog from Aurli  
# HOW TO RUN: sbatch nameofscript.sh
# DEPENDENCIES: InterproScan

###MODULES
module load blast+/2.10.0
module load interproscan/5.55.88
module load libsvm/3.25 # loads conda environment with python, libsvm and openjdk which is needed for interpro

#ENVIRONMENTS


# WORKING DIRECTORY. Should be laby/proteomics
pwd


# VARIABLE. 
python /gpfs/projects/RestGroup/agilgomez/tools/MultiLoc2/MultiLoc2/src/multiloc2_prediction.py -fasta=input_fasta/Aurli1_GeneCatalog_proteins_20120618.aa.fasta -origin=fungal -result=input_anno/multiloc_anno.txt
