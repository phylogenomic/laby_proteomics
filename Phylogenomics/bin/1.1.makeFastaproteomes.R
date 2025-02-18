library(pacman)
pacman::p_load(tidyverse,Biostrings)

# Make proteomes for each species

seq_path <- "/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/phylogenomics/data/rfdb3_2020_Laby_andAurantio.fasta"

fastafile <- readAAStringSet(seq_path)
qseqid <- names(fastafile)
sequence <- paste(fastafile)
df_seq <- data.frame(qseqid, sequence) |>
  as_tibble() |>
 separate(qseqid, into = c("protein_id", "taxonomy"), sep = "_", extra = "merge")|>
 mutate(species_id=gsub("_",".",taxonomy))
 
species_s <- table(df_seq$species_id)[table(df_seq$species_id)>10000]
sp_n <- names(species_s)

df_seq1 <- df_seq |> filter(species_id %in% sp_n) |>
    mutate(seq_name=paste0(species_id,"_",protein_id))

df_seq2 <- list()
for (i in 1:length(sp_n)){
    df_seq2[[i]] <- df_seq1 |> filter(species_id==sp_n[i])
}

# Stramenopiles
seq_path <- "/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/phylogenomics/data/rfdb3_2020_str_not_laby.fasta"
fastafile <- readAAStringSet(seq_path)
qseqid <- names(fastafile)
sequence <- paste(fastafile)

df_seq_st <- data.frame(qseqid, sequence) |>
  as_tibble() |>
 separate(qseqid, into = c("protein_id", "taxonomy"), sep = "_", extra = "merge")|>
 mutate(species_id=gsub("_",".",taxonomy))|>
 mutate(species_id= gsub("\\.$", "", species_id))|>
 mutate(species_id= gsub("sp$", "sp.$", species_id))
 
species_st <- table(df_seq_st$species_id)[table(df_seq_st$species_id)>10000]
sp_st <- names(species_st)

df_seq1_st <- df_seq_st |> filter(species_id %in% sp_st) |>
    mutate(seq_name=paste0(species_id,"_",protein_id))

df_seq2_st <- list()
for (i in 1:length(sp_n_st)){
    df_seq2_st[[i]] <- df_seq1_st |> filter(species_id==sp_st[i])
}

# Non-stramenopiles:
seq_path <- "/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/phylogenomics/data/rfdb3_2020_nonstr_only.fasta"
fastafile <- readAAStringSet(seq_path)
qseqid <- names(fastafile)
sequence <- paste(fastafile)

df_seq_nst <- data.frame(qseqid, sequence) |>
  as_tibble() |>
 separate(qseqid, into = c("protein_id", "taxonomy"), sep = "_", extra = "merge")

df_seq_nst$taxonomy |> unique() |> str_split(pattern="_") |> map_chr(~ .x[2]) |> unique()

# Select TSAR and Haptista
filt_groups <- c("Ciliophora", "Dinophyceae", "Chromerida", 
                "Apicomplexa", "Pelagophyceae", "Cercozoa", "Foraminifera","Haptista")
filt_pattern <- paste(filt_groups, collapse = "|")

df_seq_nst1 <- df_seq_nst |>
  mutate(FILT = ifelse(str_detect(taxonomy, filt_pattern), TRUE, FALSE))|>
  filter(FILT==TRUE)

table(df_seq_nst1$taxonomy)

species_nst <- table(df_seq_nst1$taxonomy)[table(df_seq_nst1$taxonomy)>10000 & table(df_seq_nst1$taxonomy)<30000]

df_seq_nst2 <- df_seq_nst1 |>
mutate(species_id=taxonomy) |>
 filter(species_id %in% names(species_nst)) |>
 mutate(species_id=gsub("_",".",taxonomy))|>
 mutate(species_id= gsub("\\.$", "", species_id))|>
 mutate(species_id= gsub("sp$", "sp.$", species_id)) |> 
 mutate(seq_name=paste0(species_id,"_",protein_id)) 


df_seq2_nst <- list()
for (i in 1:length(unique(df_seq_nst2$species_id))){
 sp_name <- unique(df_seq_nst2$species_id)[i]
 df_seq2_nst[[i]] <- df_seq_nst2 |> 
    filter(species_id==sp_name)
}

# Combine sequences
df_00 <- df_seq2
df <- list()
for (j in seq_along(df_00)) {

 fileName <- df_00[[j]]$species_id |> unique()
 print(paste(j,"of",length(df_00)))

 df[[j]] <- df_00[[j]] |> select(seq_name,sequence)|>
 mutate(seq_name=paste0(">",seq_name))

  d <- do.call(rbind, lapply(seq(nrow(df[[j]])),
                             function(i) t(df[[j]][i, ])))
  write.table(d,
              file.path("/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/phylogenomics/data/proteomes",
              paste0(fileName,".fasta")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

df_00 <- df_seq2_st
df <- list()
for (j in seq_along(df_00)) {

 fileName <- df_00[[j]]$species_id |> unique()
 print(paste(j,"of",length(df_00)))

 df[[j]] <- df_00[[j]] |> select(seq_name,sequence)|>
 mutate(seq_name=paste0(">",seq_name))

  d <- do.call(rbind, lapply(seq(nrow(df[[j]])),
                             function(i) t(df[[j]][i, ])))
  write.table(d,
              file.path("/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/phylogenomics/data/proteomes",
              paste0(fileName,".fasta")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}


df_00 <- df_seq2_nst
df <- list()
for (j in seq_along(df_00)) {

 fileName <- df_00[[j]]$species_id |> unique()
 print(paste(j,"of",length(df_00)))

 df[[j]] <- df_00[[j]] |> select(seq_name,sequence)|>
 mutate(seq_name=paste0(">",seq_name))

  d <- do.call(rbind, lapply(seq(nrow(df[[j]])),
                             function(i) t(df[[j]][i, ])))
  write.table(d,
              file.path("/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/phylogenomics/data/proteomes",
              paste0(fileName,".fasta")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

