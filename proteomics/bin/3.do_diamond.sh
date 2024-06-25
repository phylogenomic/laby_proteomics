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

# Example one, blast Aurli to Stramenopiles
# sbatch --export=input=input_fasta/Aurli1_aa.fasta,output=output_blasts/Aurliprot_Stram/Aurliprot_Stram.out bin/3.do_diamond.sh

# Example two, detected from all set to JGI.
# sbatch --export=input=DEP_results_fasta/detected.all_set.fasta,output=output_blasts/diamond_out4/detected_all_blast.out bin/3.do_diamond.sh

input_db=/gpfs/projects/CollierGroup/agilgomez/projects/laby/proteomics/input_fasta/mmetsp_uniprot.fa
db=/gpfs/projects/CollierGroup/agilgomez/projects/laby/proteomics/input_fasta/mmetsp

diamond makedb  --in $input_db \
    --db $db

input_query=/gpfs/projects/CollierGroup/agilgomez/projects/laby/proteomics/input_fasta/mmetsp_uniprot.fa
output=/gpfs/projects/CollierGroup/agilgomez/projects/laby/proteomics/output_blasts/mmetsp_to_mmetsp/mmetsp_to_mmetsp.out

diamond blastp --db $db \
--query $input_query \
--very-sensitive \
--outfmt 6 qseqid sseqid pident length mismatch gapopen evalue bitscore pident stitle \
--evalue 1e-5 \
--max-hsps 1 \
--max-target-seqs 10 \
--threads 96  \
--out $output