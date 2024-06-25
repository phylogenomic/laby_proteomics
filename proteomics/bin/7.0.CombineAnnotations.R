library(pacman)

pacman::p_load(tidyverse, Biostrings, FactoMineR,
               data.table, plotly, DEP, sva,
               SummarizedExperiment, ComplexHeatmap,
               patchwork, factoextra,
               nVennR, ggrepel, RColorBrewer, ggVennDiagram, grid, here)


# Add annotations:

# Orthoclusters
ortho <- read_delim("proteomics/input_anno/keggAu.txt", col_names = FALSE)
colnames(ortho)[1] <- "name_anno"
colnames(ortho)[7] <- "orthocluster"
ortho <- ortho |>
  select(-c("X2", "X3", "X4", "X5", "X6", "X8")) |>
  mutate(name_anno = as.character(name_anno))

#What clusters are present in at least 3 taxa?
clust_taxa <- read_csv("proteomics/input_anno/Honfer1-comparative-qc.4746.csv") |>
  mutate(cluster = as.numeric(cluster))
clust_taxa <- clust_taxa[-15367, -7]
colnames(clust_taxa) <- c("orthocluster", "Aplke1", "Aurli1",
                          "Honfer1", "LabyHa", "Schag1")

orthoclust <- ortho |>
  left_join(clust_taxa, by = "orthocluster") |>
  mutate(
      Aplke1L = Aplke1 != 0,
      Aurli1L = Aurli1 != 0,
      Honfer1L = Honfer1 != 0,
      LabyHaL = LabyHa != 0,
      Schag1L = Schag1 != 0) |>
  rowwise() |>
  mutate(num.species = sum(Aplke1L, Aurli1L, Honfer1L, LabyHaL, Schag1L),
      conserved = num.species >= 3) |>
  ungroup() |> 
  select(-ends_with("L"))

# Annotation files.
# Diamond annotations.
anno <- read_csv("proteomics/output_blasts/aurli_to_aurli/jgiblastannotations.csv") |>
  rename("name_anno" = "JGI_ID") |>
  mutate(name_anno = as.character(name_anno)) |>
  group_by(name_anno) |>
  summarise(Annotations = paste0(Annotations, collapse = "; ")) |> 
  mutate(genbank=Annotations) |> 
  separate(genbank,into=LETTERS[1:5],sep=" ") |> 
  select(name_anno,Annotations,A) |> 
  rename(genbank=A)

anno0 <- read_delim("proteomics/input_anno/Hondaea_uniprot_genbank.tsv") |> 
  select(From,Entry) |> 
  rename(Uniprot_Hon=Entry)|> 
  rename(genbank=From)

anno <- left_join(anno,anno0,by="genbank")


# GO terms
anno1 <- read_tsv(
  "GO_comparison/Aurli1_GeneCatalog_proteins_20120618_GO.tab") |>
  rename(name_anno = "#proteinId") |>
  mutate(name_anno = as.character(name_anno)) |>
  group_by(name_anno) |>
  summarise(gotermId = paste0(gotermId, collapse = "; "),
            goName = paste0(goName, collapse = "; "),
            gotermType = paste0(gotermType, collapse = "; "),
            goAcc = paste0(goAcc, collapse = "; "))
# kogdefline, kogClass
anno2 <- read_delim(
  "proteomics/input_anno/aurli1-genecatalog-proteins-20120618-kog.tab") |>
  select(proteinId, kogdefline, kogClass) |>
  rename(name_anno = proteinId) |>
  mutate(name_anno = as.character(name_anno)) |>
  mutate(kogClass = gsub("[[:space:]]*$", "", kogClass)) |>
  mutate(kogClass = ifelse(kogClass == "no description",
                           "Function unknown", kogClass)) |>
  group_by(name_anno) |>
  summarise(kogdefline = paste0(kogdefline, collapse = "; "),
            kogClass = paste0(kogClass, collapse = "; "))
# GO compartment
anno3 <- read_tsv("GO_comparison/GO_predi_det_Aurli1_ekhidna2.txt") |>
  separate(qpid, into = c("a", "b", "c", "d")) |>
  rename(name_anno = c) |>
  select(-c(a, b, d)) |>
  select(name_anno, desc) |>
  group_by(name_anno) |>
  summarise(GOcompartment = paste0(desc, collapse = "; "))
# Significant Pfams
anno4 <- read_csv("proteomics/input_anno/Aurli_sigpfams.csv") |>  
  select(2:4) |> 
  rename(name_anno = V1,
         Sig_Pfam_info = V2,
         Sig_Pfam_number = V3) |>
  separate(name_anno, sep = "\\|", into = letters[1:4]) |>
  select(-c("a", "b", "d")) |>
  rename("name_anno" = "c") |>
  mutate(name_anno = as.character(name_anno)) |>
  group_by(name_anno) |>
  summarise(Sig_Pfam_info = paste0(Sig_Pfam_info, collapse = "; "),
            Sig_Pfam_number = paste0(Sig_Pfam_number, collapse = "; "))
# InterproScan annotations.
anno6 <- read_tsv(
  "proteomics/input_anno/interpro_output.txt.tsv",col_names = FALSE) |> 
  separate(X1, sep = "\\|", into = letters[1:4]) |>
  select(-c("a", "b", "d")) |>
  rename("name_anno" = "c") |>
  mutate(name_anno = as.character(name_anno)) |>
  group_by(name_anno) |> 
  select(-X2) |> 
  group_by(name_anno,X4) |> 
  summarise("X5_X6" = paste(X5,X6, collapse = "; "),
            "X12_X13" = paste(X12,X13, collapse = "; ")) |> 
  rename(type=X4,info=X5_X6,InterproScan=X12_X13) 

anno6.1 <- anno6 |> 
  filter(type=="Pfam") |> 
  rename(Pfam_info=info,Pfam_IPS=InterproScan) |> 
  select(-type)
anno6.2 <- anno6 |> 
  filter(type=="PANTHER")|> 
  rename(PANTHER_info=info,PANTHER_IPS=InterproScan) |> 
  select(-type)

# MultiLoc2.
anno7 <- read_tsv(
  "proteomics/input_anno/multiloc_anno.txt",col_names = FALSE) |> 
  separate(X1, sep = "\\|", into = letters[1:4]) |>
  select(-c("a", "b", "d")) |>
  rename("name_anno" = "c") |>
  mutate(name_anno = as.character(name_anno)) |> 
  select(name_anno,X2,X3,X4) |> 
  rename(multiloc1=X2,multiloc2=X3,multiloc3=X4)
#Kog WGA. JGI KOG annotation.
anno8 <- read_tsv(
  "proteomics/input_anno/kog/kog-raw.txt",col_names = FALSE)|> 
  separate(X1, sep = "\\|", into = letters[1:4]) |>
  select(-c("a", "b", "d")) |>
  rename("name_anno" = "c") |> 
  select(name_anno,X13,X15) |> 
  group_by(name_anno) |> 
  rename(KOG_number=X13,KOG_info=X15) |> 
  distinct()

#KogClasses related to actin,actin-dependent.
anno8.1 <- read_csv("proteomics/input_anno/actin_kogs.csv") |> 
  rename(KOG_number=`KOG ID`)

anno8 <- anno8 |> 
  mutate(KOG_actin=ifelse(KOG_number %in% anno8.1$KOG_number,"TRUE","FALSE"))

anno8 <- anno8|> 
  summarise("KOG_number" = paste0(KOG_number, collapse = "; "),
            "KOG_info" = paste0(KOG_info, collapse = "; "),
            "KOG_with_actin" = paste0(KOG_actin, collapse = "; "))


#Keg
anno9 <- read_tsv(
  "proteomics/input_anno/Keg_numbers.txt",col_names = FALSE) |> 
  separate(X1, sep = "\\|", into = letters[1:4]) |>
  select(-c("a", "b", "d")) |>
  rename("name_anno" = "c") |> 
  rename(Keg=X2)
#KEG_brite
anno10 <- read_tsv("proteomics/input_anno/proteome_KEGG_mapping.tsv") |> 
  select(1:3,starts_with("ko")) |> 
  pivot_longer(names_to="ko",values_to="ko_values",4:45) |> 
  select(!ko)  |> 
  filter(!is.na(ko_values)) |> 
  group_by(name) |> 
  summarise("ko_values"= paste0(ko_values, collapse = "; ")) |> 
  rename(name_anno=name)
#Unique to Aurli, Labys or Different in Stramenopiles vs Labys.
anno11 <- read_csv("proteomics/input_anno/Aurliprot_conserved_Stram.csv") |> 
  select(-1) |> 
  separate(qseqid, sep = "\\|", into = letters[1:4]) |>
  select(-c("a", "b", "d")) |>
  rename("name_anno" = "c",
         "conservation_group"="Group")




fastafile <- readAAStringSet("proteomics/input_fasta/Aurli1_aa.fa")
seq_name <- names(fastafile)
sequence <- paste(fastafile)
df <- data.frame(seq_name, sequence) |>
  as_tibble() |>
  distinct() |>
  mutate(seq_name1 = seq_name) |>
  mutate(seq_name = paste0(">", seq_name)) |>
  separate(seq_name1, into = c("a", "b", "name_anno", "d"), sep = "\\|") |>
  select(seq_name, name_anno,sequence)

# Combine annotations
results <- orthoclust |>
  left_join(anno, by = "name_anno") |>
  left_join(anno1, by = "name_anno") |>
  left_join(anno2, by = "name_anno") |>
  left_join(anno3, by = "name_anno") |>
  left_join(anno4, by = "name_anno") |>
  left_join(anno6.1, by = "name_anno") |> 
  left_join(anno6.2, by = "name_anno") |> 
  left_join(anno7, by = "name_anno") |> 
  left_join(anno8, by = "name_anno") |> 
  left_join(anno9, by = "name_anno") |> 
  left_join(anno10, by = "name_anno") |> 
  left_join(anno11, by = "name_anno") |> 
  left_join(df,by= "name_anno")

write_csv(results,"proteomics/input_anno/all_anno_combined.csv")
