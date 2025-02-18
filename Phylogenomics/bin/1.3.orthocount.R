library(pacman)
pacman::p_load(tidyverse,Biostrings,furrr)

setwd("/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/")
plan(multicore, workers = parallel::detectCores() - 1)

# Count proteins for each orthogroup

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

# Load fasta files for each orthogroup
seq_path <- "/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/phylogenomics/data/proteomes/OrthoFinder/Results_Nov11/Orthogroup_Sequences/"
files <- list.files(seq_path)
possible_files <- paste0(candidate_orthogroups,".fa")
possible_files[possible_files %in% files] |> length()
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

Aurli_df4.aurli <- Aurli_df4 %>% group_by(Orthogroup) %>% 
 summarize(Aurli_JGI = sum(str_detect(protein_id, "^Aurli")))

Aurli_df4.campe <- Aurli_df4 %>% group_by(Orthogroup) %>% 
 summarize(Aurli_CAMPEP = sum(str_detect(protein_id, "^CAMPEP")))

Aurli_df4.other <- Aurli_df4 %>% group_by(Orthogroup) %>% 
summarize(Aurli_other = sum(!str_detect(protein_id, "^CAMPEP") & !str_detect(protein_id, "^Aurli")))

Aurli_df5 <- left_join(Aurli_df4.aurli,Aurli_df4.campe,by="Orthogroup")|>
    left_join(Aurli_df4.other,by="Orthogroup")

ortho2 <- left_join(ortho1,Aurli_df5,by="Orthogroup")|>
mutate(Aurli_JGI=ifelse(is.na(Aurli_JGI),0,Aurli_JGI))|>
mutate(Aurli_CAMPEP=ifelse(is.na(Aurli_CAMPEP),0,Aurli_CAMPEP))|>
mutate(Aurli_other=ifelse(is.na(Aurli_other),0,Aurli_other))

ortho2 |>
 write_csv(file="/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/phylogenomics/data/ortholabyanno.csv")

# Save all Aurli CAMPEP as diamond database:
Aurli_CAMPEP <- df_seq |> 
  future_map(~ filter(.x, str_detect(seq_names, "Aurantio") & str_detect(seq_names, "CAMPEP")))
  
df_AurliCAMPEP <- bind_rows(Aurli_CAMPEP)

df_AurliCAMPEP <- df_AurliCAMPEP |>
 mutate(seq_names=paste0(">",seq_names))

d <- do.call(rbind, lapply(seq(nrow(df_AurliCAMPEP)),
                           function(i) t(df_AurliCAMPEP[i, ])))
write.table(d,
              file.path("/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/phylogenomics/data/",
              "AurliCAMPEP.fasta"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
