# Blast the query mmetsp_uniprot.fa, to db of all Aurli1_GeneCatalog_proteins.

formatdb -t Aurlijgi -i input_fasta/Aurli1_GeneCatalog_proteins_20120618.aa.fasta
blastp -query input_fasta/mmetsp_uniprot.fa \
  -db input_fasta/Aurli1_GeneCatalog_proteins_20120618.aa.fasta \
  -out input_anno/blast_mmetsptojgi.out -outfmt 6 
  
# There is no E.value filter, it is filter later for 10^20

  
#to check progress:
awk -F '|' '{print $1}'  mmetsp-uniprot-to-jgi.out | sort | uniq | wc -l
