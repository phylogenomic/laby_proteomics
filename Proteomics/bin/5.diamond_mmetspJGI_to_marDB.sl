#!/usr/bin/env bash
#SBATCH --job-name=diamond_nr_long96
#SBATCH --output=diamond_nr.log
#SBATCH --ntasks-per-node=96
#SBATCH --nodes=1
#SBATCH --time=2-00:00:00
#SBATCH -p long-96core


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

# TEST
#sbatch --export=input_db='/gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020.fasta',
# input_query='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/input_fasta/query.fa',
# output_blast='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/output_blasts/diamond_to_marDB/test.out' bin/5.diamond_mmetspJGI_to_marDB.sl

module load hts/1.0

# DB0. Aurli only
seqkit grep -n -r -p "Aurantiochytrium_limacinum" /gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020.fasta > /gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020_Aurli_only.fasta
input_db='/gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020_Aurli_only.fasta'

# DB1. Aurantio without Aurli
seqkit grep -n -r -f /gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/output_blasts/diamond_to_marDB/Aurantio_notAurli.txt /gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020.fasta > /gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020_Aurantio_notAurli.fasta
input_db='/gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020_Aurantio_notAurli.fasta'

# DB2. Labys without Aurantio
seqkit grep -n -r -f /gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/output_blasts/diamond_to_marDB/Laby_notAurantio.txt /gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020.fasta > /gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020_Laby_notAurantio.fasta
input_db='/gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020_Laby_notAurantio.fasta'

# DB3. Non-stramenopiles
seqkit grep -n -r -f /gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/output_blasts/diamond_to_marDB/nonstramenopile_only.txt /gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020.fasta -o /gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020_nonstr_only.fasta
input_db='/gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020_nonstr_only.fasta'

# DB4. Stramenopiles without Labys or Aurantio
seqkit grep -n -r -f /gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/output_blasts/diamond_to_marDB/stramenopile_only.txt /gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020.fasta |
seqkit grep -n -r -v -f /gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/output_blasts/diamond_to_marDB/Aurantio_and_other_laby.txt -o /gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020_str_not_laby.fasta
input_db='/gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020_str_not_laby.fasta'

# Make DB
module load diamond/2.0.10
input_db='/gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020_Aurli_only.fasta'
input_db='/gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020_Aurantio_notAurli.fasta'
input_db='/gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020_Laby_notAurantio.fasta'
input_db='/gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020_nonstr_only.fasta'
input_db='/gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020_str_not_laby.fasta'
diamond makedb  --in $input_db \
    --db $input_db.db


# Search:
# DB0
input_db='/gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020_Aurli_only.fasta'
input_query='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/input_fasta/ox87102.2020_06.faa'
output_blast='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/output_blasts/diamond_to_marDB/mmetspENA_to_Aurli_only.out'

input_query='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/input_fasta/mmetsp_uniprot.fa'
output_blast='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/output_blasts/diamond_to_marDB/mmetsp_to_Aurli_only.out'

input_query='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/input_fasta/Aurli1_GeneCatalog_proteins_20120618.aa.fasta'
output_blast='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/output_blasts/diamond_to_marDB/JGI_to_Aurli_only.out'

diamond blastp --db $input_db.db \
--query $input_query \
--very-sensitive \
--outfmt 6 qseqid sseqid pident length mismatch gapopen evalue bitscore stitle \
--evalue 1e-3 \
--max-hsps 1 \
--max-target-seqs 50 \
--threads 96  \
--out $output_blast # \ --masking 0 No-masking


#DB1
input_db='/gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020_Aurantio_notAurli.fasta'
input_query='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/input_fasta/ox87102.2020_06.faa'
output_blast='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/output_blasts/diamond_to_marDB/mmetspENA_to_Aurantio_notAurli.out'

input_query='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/input_fasta/mmetsp_uniprot.fa'
output_blast='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/output_blasts/diamond_to_marDB/mmetsp_to_Aurantio_notAurli.out'

input_query='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/input_fasta/Aurli1_GeneCatalog_proteins_20120618.aa.fasta'
output_blast='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/output_blasts/diamond_to_marDB/JGI_to_Aurantio_notAurli.out'

diamond blastp --db $input_db.db \
--query $input_query \
--very-sensitive \
--outfmt 6 qseqid sseqid pident length mismatch gapopen evalue bitscore stitle \
--evalue 1e-3 \
--max-hsps 1 \
--max-target-seqs 50 \
--threads 96  \
--out $output_blast # \ --masking 0 No-masking

# DB2
input_db='/gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020_Laby_notAurantio.fasta'
input_query='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/input_fasta/ox87102.2020_06.faa'
output_blast='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/output_blasts/diamond_to_marDB/mmetspENA_to_Laby_notAurantio.out'

input_query='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/input_fasta/mmetsp_uniprot.fa'
output_blast='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/output_blasts/diamond_to_marDB/mmetsp_to_Laby_notAurantio.out'

input_query='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/input_fasta/Aurli1_GeneCatalog_proteins_20120618.aa.fasta'
output_blast='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/output_blasts/diamond_to_marDB/JGI_to_Laby_notAurantio.out'

diamond blastp --db $input_db.db \
--query $input_query \
--very-sensitive \
--outfmt 6 qseqid sseqid pident length mismatch gapopen evalue bitscore stitle \
--evalue 1e-3 \
--max-hsps 1 \
--max-target-seqs 10 \
--threads 96  \
--out $output_blast # \ --masking 0 No-masking

# DB3
input_db='/gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020_nonstr_only.fasta'
input_query='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/input_fasta/ox87102.2020_06.faa'
output_blast='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/output_blasts/diamond_to_marDB/mmetspENA_to_nonstr_only.out'

input_query='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/input_fasta/mmetsp_uniprot.fa'
output_blast='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/output_blasts/diamond_to_marDB/mmetsp_to_nonstr_only.out'

input_query='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/input_fasta/Aurli1_GeneCatalog_proteins_20120618.aa.fasta'
output_blast='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/output_blasts/diamond_to_marDB/JGI_to_nonstr_only.out'

diamond blastp --db $input_db.db \
--query $input_query \
--very-sensitive \
--outfmt 6 qseqid sseqid pident length mismatch gapopen evalue bitscore stitle \
--evalue 1e-3 \
--max-hsps 1 \
--max-target-seqs 10 \
--threads 96  \
--out $output_blast # \ --masking 0 No-masking

# DB4
input_db='/gpfs/projects/RestGroup/mariana/carot/ref/rfdb3_2020_str_not_laby.fasta'
input_query='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/input_fasta/ox87102.2020_06.faa'
output_blast='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/output_blasts/diamond_to_marDB/mmetspENA_to_str_not_laby.out'

input_query='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/input_fasta/mmetsp_uniprot.fa'
output_blast='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/output_blasts/diamond_to_marDB/mmetsp_to_str_not_laby.out'

input_query='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/input_fasta/Aurli1_GeneCatalog_proteins_20120618.aa.fasta'
output_blast='/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/proteomics/output_blasts/diamond_to_marDB/JGI_to_str_not_laby.out'

diamond blastp --db $input_db.db \
--query $input_query \
--very-sensitive \
--outfmt 6 qseqid sseqid pident length mismatch gapopen evalue bitscore stitle \
--evalue 1e-3 \
--max-hsps 1 \
--max-target-seqs 10 \
--threads 96  \
--out $output_blast # \ --masking 0 No-masking