#!/usr/bin/env Rscript
library(pacman)

pacman::p_load(tidyverse, data.table,plotly,here,grid,patchwork,ggrepel,
                RColorBrewer,ggVennDiagram,factoextra,FactoMineR,furrr)# CRAN
pacman::p_load(Biostrings)

plan(multicore, workers = parallel::detectCores() - 1)
setwd("/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/")
options(future.globals.maxSize = 800 * 1024^2)

# Orthofinder results
ortho <- read.table("phylogenomics/data/proteomes/OrthoFinder/Results_Nov11/Orthogroups/Orthogroups.txt",
header=FALSE,sep=":")
 
ortho <- ortho |>
  dplyr::rename(orthogroup=V1) |>
  dplyr::rename(allspecies=V2)

results <- read_csv("proteomics/input_anno/all_anno_combined.csv")


prot_ortho <- future_map(1:nrow(results), function(i) {
  prot_num <- results$name_anno[i]
  prot1 <- paste0("Aurantiochytrium.limacinum_Aurli", prot_num," ")
  prot2 <- paste0("Aurantiochytrium.limacinum_Aurli", prot_num,"$")
  print(paste0(i, " of ", nrow(results)))

  orthogroup <- ortho %>%
    filter(str_detect(allspecies, prot1)|str_detect(allspecies, prot2)) %>%
    select(orthogroup) %>%
    pull()

  result <- paste0(prot_num, "-", orthogroup)
  print(result)
  return(result)
})

prot_ortho <- unlist(prot_ortho)

split_elements <- strsplit(prot_ortho, "-")
df <- data.frame(
   name_anno = sapply(split_elements, function(x) x[1]),
   orthogroup = sapply(split_elements, function(x) x[2])
 ) 
 
df <- df |>
  distinct()

df |> filter(is.na(orthogroup)) |> dim() # 0
df |> dim() # 14859

# qseqid and Orthogroup
df <- df |>
mutate(Orthogroup=orthogroup,
        qseqid=name_anno)|>
        select(qseqid,Orthogroup)

write_csv(df,file="phylogenomics/data/jgi_to_orthogroup.csv")
