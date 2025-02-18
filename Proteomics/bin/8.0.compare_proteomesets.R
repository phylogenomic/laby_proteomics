library(pacman)

pacman::p_load(tidyverse, data.table,plotly,here,grid,patchwork,ggrepel, # CRAN
               RColorBrewer,ggVennDiagram,factoextra,FactoMineR, #CRAN
               ggplotify,ggvenn, # CRAN
               Biostrings,DEP, sva,SummarizedExperiment, ComplexHeatmap)

x <- read_delim("proteomics/output_blasts/mmetsp_and_jgi/jgi_to_mmetsp.out",col_names = FALSE)
blast_jgi2mmetsp <- x |> 
  mutate(X3=as.numeric(X3),X4=as.numeric(X4),
         X5=as.numeric(X5),X6=as.numeric(X6),
         X7=as.numeric(X7),X8=as.numeric(X8),
         X9=as.numeric(X9)) |> 
  group_by(X1) |> 
  slice_min(n = 1, X7) |> #E.value
  slice_max(n = 1, X8) |> #bit
  slice_max(n = 1, X4) |> #length
  slice_max(n = 1, X3) |> #pident
  slice_max(n = 1, X9) |> #pident
  mutate(X1=gsub(';','',X1)) |> 
  slice_sample(n=1)|> 
  rename(seq_name1=X1)|> 
  separate(seq_name1, into = c("a", "b", "name_anno", "d"), sep = "\\|") |> 
  ungroup() |> 
  mutate(ID=paste0(X2,":",name_anno)) |> 
  select(ID,X7,X8,X4,X9)|> 
  mutate(col="blast_jgi2mmetsp")
  
y <- read_delim("proteomics/output_blasts/mmetsp_and_jgi/mmetsp_to_jgi.out",col_names = FALSE)
blast_mmetsp2jgi <- y |> 
  mutate(X3=as.numeric(X3),X4=as.numeric(X4),
         X5=as.numeric(X5),X6=as.numeric(X6),
         X7=as.numeric(X7),X8=as.numeric(X8),
         X9=as.numeric(X9)) |> 
  group_by(X1) |> 
  slice_min(n = 1, X7) |> #E.value
  slice_max(n = 1, X8) |> #bit
  slice_max(n = 1, X4) |> #length
  slice_max(n = 1, X3) |> #pident
  slice_max(n = 1, X9) |> #pident
  mutate(X1=gsub(';','',X1)) |> 
  slice_sample(n=1) |> 
  rename(seq_name1=X2) |> 
  separate(seq_name1, into = c("a", "b", "name_anno", "d"), sep = "\\|")|> 
  ungroup() |> 
  mutate(ID=paste0(X1,":",name_anno)) |> 
  select(ID,X7,X8,X4,X9) |> 
  mutate(col="blast_mmetsp2jgi")

# Create non-redundant list. How many proteins are there in the non-redundant list?
# How many proteins are missing in the JGI/MMETSP?

merged_set <- full_join(blast_mmetsp2jgi,blast_jgi2mmetsp,by="ID")

# JGI has 14895 proteins
# MMETSP has 14622 proteins
14859-11503 #JGI no hits. 3356
14622-12736 #MMETSP no hits. 1886

# 14602+3356+1886=19844

# 14602 after merged (full_join).

# 9637 (inner_join). Total bidirectional best hits.

# What is the total non-redundant proteome size?

# filter by e.value




