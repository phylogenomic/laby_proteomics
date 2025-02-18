#!/bin/sh

#SBATCH --job-name=interpros
#SBATCH -o interpro.sh.out
#SBATCH --nodes=1
#SBATCH -p extended-40core
#SBATCH --time=7-00:00:00

# AUTHOR: ALEX GIL GOMEZ
# INPUTS: GeneCatalog
# OUTPUT: Pfam
# INFO: Runs interproscan on gene catalog from Aurli  
# HOW TO RUN: sbatch nameofscript.sh
# DEPENDENCIES: InterproScan

###MODULES
module load blast+/2.10.0
module load interproscan/5.55.88
module load anaconda/3
source activate openjdk
#ENVIRONMENTS


# WORKING DIRECTORY. Should be laby/proteomics
pwd


# VARIABLE. 
interproscan.sh -i input_fasta/Aurli1_GeneCatalog_proteins_20120618.aa.fasta -f tsv -b input_anno/interpro_output.txt
