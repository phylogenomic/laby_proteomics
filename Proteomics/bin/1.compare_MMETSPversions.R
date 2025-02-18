library(pacman)

pacman::p_load(tidyverse, data.table,plotly,here,grid,patchwork,ggrepel, # CRAN
               RColorBrewer,ggVennDiagram,factoextra,FactoMineR, #CRAN
               Biostrings,DEP, sva,SummarizedExperiment, ComplexHeatmap)

# Comparison of different MMETSP versions, to check what Creative Proteomics used.


fastafile <- readAAStringSet("proteomics/input_fasta/Aurli1_GeneCatalog_proteins_20120618.aa.fa",
                             format="fasta")
seq_name <- names(fastafile)
sequence <- paste(fastafile)
df1 <- data.frame(seq_name, sequence) |>
  as_tibble() |>
  distinct() |>
  mutate(seq_name = paste0(">", seq_name)) 


fastafile <- readAAStringSet("proteomics/input_fasta/mmetsp_uniprot.fa.fa")
seq_name <- names(fastafile)
sequence <- paste(fastafile)
df2 <- data.frame(seq_name, sequence) |>
  as_tibble() |>
  distinct() |>
  mutate(seq_name = paste0(">", seq_name))


a <- left_join(df2,df1,by="sequence")
a |> select(seq_name.x,seq_name.y,sequence) |> 
  write.csv("proteomics/input_anno/mmetspseqs.csv",row.names = FALSE)
