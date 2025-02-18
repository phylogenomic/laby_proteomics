library(pacman)
pacman::p_load(tidyverse,Biostrings,furrr)
plan(multicore, workers = parallel::detectCores() - 1)

# This script takes as input the MMETSP proteome to different AurliCAMPEP from Mariana's DB.

setwd("/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/")

diam_colnames <- c("qseqid","sseqid","pident","length","mismatch","gapopen",
                "evalue","bitscore","stitle")
CAMPEP_only <- read_table("phylogenomics/data/mmetsp_to_CAMPEP.out",
                    col_names = FALSE)
colnames(CAMPEP_only) <- diam_colnames

best_diamond <- CAMPEP_only |>
group_by(qseqid)|>
  slice_min(order_by = evalue,n = 1) |> 
  slice_max(order_by = pident,n = 1)|> 
  slice_max(order_by = bitscore,n = 1)|> 
  slice_max(order_by = length,n = 1)|> 
  slice_min(order_by = mismatch,n = 1)|> 
  slice_min(order_by = gapopen,n = 1)|>
  separate(stitle,into=c("taxonomy","protein_id"),sep="_")

# Get orthogroups:
ortho <- read.delim("phylogenomics/data/proteomes/OrthoFinder/Results_Nov11/Orthogroups/Orthogroups.GeneCount.tsv")|>
  mutate(across(starts_with("Eukaryota"), as.numeric))

ortho1 <- ortho %>% mutate( Total_Laby_genes = rowSums(ortho[, 66:73], na.rm = TRUE), 
Total_NotLaby_genes = rowSums(ortho[, c(2:65, 74:142)], na.rm = TRUE), 
Total_Oomycetes_genes = rowSums(ortho[, 126:133], na.rm = TRUE), 
Total_Laby_sp = rowSums(ortho[, 66:73] > 1, na.rm = TRUE), 
Total_notLaby_sp = rowSums(ortho[, c(2:65, 74:142)] > 1, na.rm = TRUE), 
Total_Oomycetes_sp = rowSums(ortho[, 126:133] > 1, na.rm = TRUE) ) |>
    mutate(Perc_Laby_genes=round(100*Total_Laby_genes/Total,2))|>
    mutate(Perc_notLaby_genes=round(100*Total_NotLaby_genes/Total,2))|>
    mutate(Perc_Oomycetes_genes=round(100*Total_Oomycetes_genes/Total,2)) |>
    arrange(Total_Laby_sp)

ortho1$Total |> sum() # 4849469
ortho1$Total_Laby_sp |> sum() # 48021
ortho1 |>
select(Eukaryota.Bigyra.Labyrinthulomycetes.Thraustochytriaceae.Aurantiochytrium.limacinum) |> 
sum() # 65838

candidate_orthogroups <- unique(ortho1$Orthogroup)

seq_path <- "/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/phylogenomics/data/proteomes/OrthoFinder/Results_Nov11/Orthogroup_Sequences/"
files <- list.files(seq_path)
possible_files <- paste0(candidate_orthogroups,".fa")
possible_files[possible_files %in% files] |> length() #302461
possible_files <- paste0(seq_path,possible_files)

fastafile <- future_map(possible_files, function(path) { 
  readAAStringSet(path) })

qseqid <- map(fastafile,names)
sequence <- future_map(fastafile,paste)
df_seq <- future_map2(qseqid, sequence,data.frame) |>
  future_map(as_tibble)

df_seq <- map(df_seq, ~ setNames(.x, c("seq_names", "seqs")))

Aurli_df <- df_seq |> 
  future_map(~ filter(.x, str_detect(seq_names, "Aurantio")))
Aurli_df2 <- future_map2(Aurli_df, candidate_orthogroups, ~ mutate(.x, Orthogroup = .y))
Aurli_df3 <- Aurli_df2 |>
  future_map(~ separate(.x, seq_names, into = c("taxonomy", "protein_id"), sep = "_"))
Aurli_df4 <- bind_rows(Aurli_df3)
Aurli_df4 <- Aurli_df4 |> filter(str_detect(protein_id, "^CAMPEP"))


# Combine MMTESP results with Orthogroup information:
orthoprot <- Aurli_df4 |> select(protein_id,Orthogroup)
best_diamond_ortho <- best_diamond |> left_join(orthoprot,by="protein_id")

best_diamond_ortho |> select(qseqid,protein_id,Orthogroup) |> distinct() 
best_diamond_ortho |> select(qseqid,Orthogroup) |> distinct() 

result <- best_diamond_ortho |>
 select(qseqid,Orthogroup)|>distinct()|>
  group_by(qseqid) |>
  summarise(Orthogroup = paste(Orthogroup, collapse = ", ")) #3252

write_csv(result,file="phylogenomics/data/mmetsp_to_orthogroup.csv")
