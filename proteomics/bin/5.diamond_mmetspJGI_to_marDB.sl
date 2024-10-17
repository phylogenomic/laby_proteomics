#!/usr/bin/env bash
#SBATCH --job-name=diamond_nr_long96
#SBATCH --output=diamond_nr.log
#SBATCH --ntasks-per-node=96
#SBATCH --nodes=1
#SBATCH --time=2-00:00:00
#SBATCH -p long-96core

module load diamond/2.0.10

# diamond blast proteins against previously formated nr db to get annotations

# running in the 40-core queues
# using "--very-sensitive" which is the second most sensitive option
# see https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#sensitivity-modes for sensitivity options


#RUN WITH:

# MMETSP (UNIPARC)
#sbatch --export=input_db='/gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020.fasta',
# input_query='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/input_fasta/mmetsp_uniprot.fa',
# output_blast='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/output_blasts/diamond_to_marDB/mmetsp_to_marDB.out' bin/5.diamond_mmetspJGI_to_marDB.sl

# MMETSP (ENA)
#sbatch --export=input_db='/gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020.fasta',
# input_query='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/input_fasta/ox87102.2020_06.faa',
# output_blast='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/output_blasts/diamond_to_marDB/mmetsp_ENA_to_marDB.out' bin/5.diamond_mmetspJGI_to_marDB.sl

# JGI
#sbatch --export=input_db='/gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020.fasta',
# input_query='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/input_fasta/Aurli1_GeneCatalog_proteins_20120618.aa.fasta',
# output_blast='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/output_blasts/diamond_to_marDB/JGI_to_marDB.out' bin/5.diamond_mmetspJGI_to_marDB.sl

# Mariana's db: /gpfs/projects/RestGroup/mariana/carot/ref
# Aurli28575_Eukaryota_Bigyra_Labyrinthulomycetes_Thraustochytriaceae_Aurantiochytrium_limacinum

diamond makedb  --in $input_db \
    --db $input_db.db

diamond blastp --db $input_db.db \
--query $input_query \
--very-sensitive \
--outfmt 6 qseqid sseqid pident length mismatch gapopen evalue bitscore stitle \
--evalue 1e-3 \
--max-hsps 1 \
--max-target-seqs 10 \
--threads 96  \
--out $output_blast