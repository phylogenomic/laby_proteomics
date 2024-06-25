library(pacman)

pacman::p_load(tidyverse,Biostrings,FactoMineR,
               data.table,plotly,DEP,sva,
               SummarizedExperiment,ComplexHeatmap,
               patchwork,factoextra,plotly,
               ggrepel,RColorBrewer,ggVennDiagram,grid,igraph)


#JGI
filename <- "rawnormalized_revised.csv"

raw_data <- read_csv(file.path("proteomics/input_data", filename))
colnames(raw_data) <- c("name_ID", "Fasta_headers",
                        "RAW.T0-1", "RAW.T0-2", "RAW.T0-3",
                        "RAW.T2-1", "RAW.T2-2", "RAW.T2-3",
                        "RAW.T4-1", "RAW.T4-2", "RAW.T4-3",
                        "RAW.T6-1", "RAW.T6-2", "RAW.T6-3",
                        "RAW.T8-1", "RAW.T8-2", "RAW.T8-3")
jgi_IDs <- raw_data |>
  separate(name_ID, sep = "\\|", into = letters[1:4]) |>
  select(-c("a", "b", "d")) |>
  rename("name_ID" = "c") |>
  select(-Fasta_headers) |> 
  select(name_ID) |> distinct() |> pull()


# Blast results JGI-MMETSP mapping:
# Earlier results
#old_blast <- "proteomics/input_anno/blast_mmetsptojgi.out"
new_blast <- "proteomics/input_anno/blast_JGI_to_MMETSP_May24_filtered_results.tsv"


# Select JGI detected experimentally and get the best match using filters.


# For MMETSP
raw_blast <- read_delim(new_blast,col_names = FALSE)

result <- raw_blast |> 
  mutate(X3=as.numeric(X3),X4=as.numeric(X4),
         X5=as.numeric(X5),X6=as.numeric(X6),
         X7=as.numeric(X7),X8=as.numeric(X8),
         X9=as.numeric(X9),X10=as.numeric(X10),
         X11=as.numeric(X11),X12=as.numeric(X12)) |> 
  group_by(X1)|> 
  separate(X2, into = c("a", "b", "JGI_ID", "d"), sep = "\\|") |> 
  select(-c("a","b","d")) |> 
  mutate(present_MQ=ifelse(JGI_ID %in% jgi_IDs,TRUE,FALSE)) |> 
  rename(MMETSP_ID=X1,Percent_Identity=X3,
       Aminoacid_length=X4,E.value=X11,Bitscore=X12)


# TRUE are JGI present in MQ
result1_TRUE <- result |> 
  group_by(MMETSP_ID,JGI_ID) |> 
  filter(present_MQ==TRUE) |> 
  slice_min(n = 1,E.value) |> #E value
  slice_max(n = 1, Bitscore) |> #Bitscore
  slice_max(n = 1, Aminoacid_length) |> #Amino acid length
  slice_max(n = 1, Percent_Identity) |> #Highest percent identity
  slice_max(n = 1, X9) |> 
  filter(E.value<10^-20)  # Filter E-values lower than 10^-20 (removes 0 proteins)
  
# FALSE are JGI not present in MQ
result1_FALSE <- result |> 
  group_by(MMETSP_ID,JGI_ID) |> 
  filter(present_MQ==FALSE) |> 
  slice_min(n = 1,E.value) |> #E value
  slice_max(n = 1, Bitscore) |> #Bitscore
  slice_max(n = 1, Aminoacid_length) |> #Amino acid length
  slice_max(n = 1, Percent_Identity) |> #Highest percent identity
  slice_max(n = 1, X9) |> 
  filter(E.value<10^-20)  # Filter E-values lower than 10^-20 (removes 0 proteins)

result1_TRUE$MMETSP_ID |> unique() |> length() #2559 MMETSP in JGI set
result1_FALSE$MMETSP_ID |> unique() |> length() #437 MMETSP not in JGI set


# If present in MQ with blast hits not present in MQ, keep best present in MQ
MQ_with_blast <- intersect(result1_FALSE$MMETSP_ID,result1_TRUE$MMETSP_ID) #1
MQ_both <- result1_TRUE |> filter(MMETSP_ID %in% MQ_with_blast)|> 
  ungroup() |> 
  group_by(MMETSP_ID) |> 
  slice_min(n = 1,E.value) |> #E value
  slice_max(n = 1, Bitscore) |> #Bitscore
  slice_max(n = 1, Aminoacid_length) |> #Amino acid length
  slice_max(n = 1, Percent_Identity) |> #Highest percent identity
  slice_max(n = 1, X9) |> 
  arrange(MMETSP_ID)|> 
  sample_n(1)

# If present in MQ without blast hits not present in MQ, keep best MQ
MQ_without_blast <- setdiff(result1_TRUE$MMETSP_ID,result1_FALSE$MMETSP_ID) 
MQ_without_blast |> unique() |> length() #2558
MQ1 <- result1_TRUE |> 
  filter(MMETSP_ID %in% MQ_without_blast)|> 
  ungroup() |> 
  group_by(MMETSP_ID) |> 
  slice_min(n = 1,E.value) |> #E value
  slice_max(n = 1, Bitscore) |> #Bitscore
  slice_max(n = 1, Aminoacid_length) |> #Amino acid length
  slice_max(n = 1, Percent_Identity) |> #Highest percent identity
  slice_max(n = 1, X9) |> 
  arrange(MMETSP_ID)|> 
  sample_n(1)

# If not present in MQ with blast hits not present in MQ, keep best blast
notMQ_with_blast <- setdiff(result1_FALSE$MMETSP_ID,result1_TRUE$MMETSP_ID)
notMQ_with_blast |> unique() |> length() #436
MQ2 <- result1_FALSE |> 
  filter(MMETSP_ID %in% notMQ_with_blast)|> 
  ungroup() |> 
  group_by(MMETSP_ID) |> 
  slice_min(n = 1,E.value) |> #E value
  slice_max(n = 1, Bitscore) |> #Bitscore
  slice_max(n = 1, Aminoacid_length) |> #Amino acid length
  slice_max(n = 1, Percent_Identity) |> #Highest percent identity
  slice_max(n = 1, X9) |> 
  arrange(MMETSP_ID) |> 
  sample_n(1)


result <- rbind(MQ_both,MQ1,MQ2)


# Which MMETSP have more than 1 hit?
# "A0A6S8EWL1" (75431,137653). 137653 is in the detected set
# "A0A6S8FCD2" (141084,138522,137655). 137655 is in the detected set
# "A0A6S8G5Y6" (148640,140387,140385,8244). 148640 is in the detected set
# 2 histones and one hypothetical protein have more than one hit.
# Select the JGI that was also in the JGI set.
not_one <- result |> 
  group_by(MMETSP_ID) |> 
  summarise(con=n()) |> 
  filter(con>1) |> 
  select(MMETSP_ID) |> pull() |> 
  str_split("\\.")|>
  map_chr(1)


result1 <- result |>
  separate("MMETSP_ID",into=c("MMETSP_ID","U2"),sep="\\.") |> 
  select(-U2) |> 
  group_by(MMETSP_ID) |> 
  slice(sample(1))  # Samples one randomly.


# After every MMETSP has a single JGI
not_one <- result1 |> 
  group_by(MMETSP_ID) |> 
  summarise(con=n()) |> 
  filter(con>1) |> 
  select(MMETSP_ID) |> pull() |> 
  str_split("\\.")|>
  map_chr(1)
not_one

# Result 1 = MMETSP mapped to JGI
result1$MMETSP_ID |> unique() |> length() #3251 unique MMETSP
#new: 2995 unique MMETSP
result1$JGI_ID |> unique() |> length() #3184 unique JGI
# new: 2951 unique JGI

result1 |> select(MMETSP_ID,JGI_ID) |> 
  write_csv("proteomics/input_anno/jgimmetsp_result_may24.csv")

# Collapse proteins that have high similar percent identity into clusters.

# JGI proteins
seq_fasta <- readAAStringSet("proteomics/input_fasta/Aurli1_aa.fa")
seq_name <- names(seq_fasta)
sequence <- paste(seq_fasta)
df <- data.frame(seq_name, sequence) |>
  as_tibble() |>
  distinct() |>
  mutate(seq_name = paste0(">", seq_name)) |>
  mutate(seq_name1 = seq_name) |>
  separate(seq_name1, into = c("a", "b", "name", "d"), sep = "\\|") |>
  select(seq_name, name, sequence)

all_JGI <- df$name

# For all by all
threshold <- 94
jgi_blast <- read_csv("proteomics/output_blasts/aurli_to_aurli/diamondAurli1AllbyAll.csv") |> 
  filter(self==0)|> 
  filter(evalue<10^-20)|> 
  separate("qseqid",into=c("qseqid1","qseqid2","qseqid"),sep="\\|") |> 
  select(-c(qseqid1,qseqid2)) |> 
  separate("sseqid",into=c("sseqid1","sseqid2","sseqid"),sep="\\|") |> 
  select(-c(sseqid1,sseqid2)) |> 
  filter(pident...3>threshold)  |> 
  mutate(in_qseqid=qseqid%in%all_JGI) |> 
  mutate(in_sseqid=sseqid%in%all_JGI) |> 
  mutate(pair=in_qseqid&in_sseqid) |> 
  filter(pair==TRUE) |> 
  select(c(1,2)) |> 
  mutate(qseqid=as.numeric(qseqid),
         sseqid=as.numeric(sseqid))

for (i in 1:dim(jgi_blast)[1]){
  qs <- jgi_blast$qseqid[i]
  ss <- jgi_blast$sseqid[i]
  pair <- sort(c(qs,ss))
  jgi_blast$pair[i] <- paste0(pair,collapse = "-")
}

jgi_blast$pair |> unique()|> length()

jgi_blast <- jgi_blast |> select(pair) |> 
  distinct() |> separate(pair,into=c("p1","p2"),sep="-") |> 
  select(p1,p2) |> 
  arrange(p1)

g <- igraph::graph_from_data_frame(jgi_blast, directed = FALSE, vertices = NULL)

plot(g)
eb <- cluster_edge_betweenness(g)
grps <- split(V(g),eb$membership)


list_clusters <- list()
for (i in seq_along(grps)){
  list_clusters[[i]] <- attributes(grps[[i]])$names |> 
    as.numeric() |> sort() |> paste0(collapse="&")
}
clusters <- list_clusters |> unlist() |> as.data.frame()

df <- clusters |> 
  rename(group=`unlist(list_clusters)`) |> 
  mutate(group_cp=group) |> 
  separate(group_cp,into=c('a','b','c','d','e','f'),sep="&") |> 
  pivot_longer(names_to = "name",values_to = "JGI_ID",2:7) |> 
  select(-name) |> 
  drop_na() |> 
  select(JGI_ID,group) |> 
  arrange(group)

# df = JGI by 94% identity JGI groups
df |> 
  rename(JGI_g=group) |> write_csv("proteomics/input_anno/jgi_highpercentidentity_may24.csv")

# MMETSP:

# Collapse proteins that have high similar percent identity into clusters.

# MMETSP proteins
seq_fasta <- readAAStringSet("proteomics/input_fasta/mmetsp_uniprot.fa")
seq_name <- names(seq_fasta)
sequence <- paste(seq_fasta)
df <- data.frame(seq_name, sequence) |>
  as_tibble() |>
  distinct() |>
  mutate(seq_name1 = seq_name) |>
  mutate(seq_name = paste0(">", seq_name)) |>
  separate(seq_name1, into = c("name1", "b"), sep = ";") |>
  separate(name1, into = c("name", "b"), sep = " \\(obsolete\\)") |> 
  select(seq_name, name, sequence)

all_MMETSP <- df$name

# For all by all
threshold <- 94
mmetsp_blast <- read_delim("proteomics/output_blasts/mmetsp_to_mmetsp/mmetsp_to_mmetsp.out",
                           col_names = FALSE) |> 
  rename(qseqid=X1,sseqid=X2,pident=X3,length=X4,
         mismatch=X5,gapopen=X6,evalue=X7,bitscore=X8,pident2=X9,
         stitle=X10) |> 
  mutate(self=ifelse(qseqid==sseqid,1,0)) |> 
  filter(self==0)|> 
  filter(evalue<10^-20)|> 
  separate(qseqid, into = c("qseqid", "b"), sep = ";") |> 
  separate(sseqid, into = c("sseqid", "c"), sep = ";") |> 
  select(-c("b","c")) |> 
  mutate(in_qseqid=qseqid%in%all_MMETSP) |> 
  mutate(in_sseqid=sseqid%in%all_MMETSP) |> 
  mutate(pair=in_qseqid&in_sseqid) |> 
  filter(pair==TRUE) |> 
  select(c(1,2)) 


for (i in 1:dim(mmetsp_blast)[1]){
  qs <- mmetsp_blast$qseqid[i]
  ss <- mmetsp_blast$sseqid[i]
  pair <- sort(c(qs,ss))
  mmetsp_blast$pair[i] <- paste0(pair,collapse = "-")
}

mmetsp_blast$pair |> unique()|> length()

mmetsp_blast <- mmetsp_blast |> select(pair) |> 
  distinct() |> separate(pair,into=c("p1","p2"),sep="-") |> 
  select(p1,p2) |> 
  arrange(p1)

g <- igraph::graph_from_data_frame(mmetsp_blast,
                                   directed = FALSE,
                                   vertices = NULL)

plot(g)
eb <- cluster_edge_betweenness(g)
grps <- split(V(g),eb$membership)


list_clusters <- list()
for (i in seq_along(grps)){
  list_clusters[[i]] <- attributes(grps[[i]])$names |> 
    sort() |> paste0(collapse="&")
}
clusters <- list_clusters |> unlist() |> as.data.frame()

df <- clusters |> 
  rename(group=`unlist(list_clusters)`) |> 
  mutate(group_cp=group) |> 
  separate(group_cp,into=c('a','b','c','d','e','f'),sep="&") |> 
  pivot_longer(names_to = "name",values_to = "MMETSP_ID",2:7) |> 
  select(-name) |> 
  drop_na() |> 
  select(MMETSP_ID,group) |> 
  arrange(group)

# df = JGI by 94% identity JGI groups
df |> 
  rename(MMETSP_g=group) |> write_csv("proteomics/input_anno/mmetsp_highpercentidentity_may24.csv")

#Blast output

# For new JGI/diamond (with annotations)
blast <- read_delim("proteomics/output_blasts/aurli_to_aurli/diamondFullAurli1proteome.out",col_names = FALSE) |> 
  mutate(X3=as.numeric(X3),X4=as.numeric(X4),
         X5=as.numeric(X5),X6=as.numeric(X6),
         X7=as.numeric(X7),X8=as.numeric(X8),
         X9=as.numeric(X9)) |> 
  group_by(X1) |> 
  slice_min(n = 1,X7) |> #E value 
  filter(X7<10^-20) |>  # Filter E-values lower than 10^-20 (removes 3304-3257= 47 proteins)
  slice_max(n = 1, X8) |> #Bitscore
  slice_max(n = 1, X4) |> #Amino acid length
  slice_max(n = 1, X3) |> #Highest percent identity
  slice_max(n = 1, X9) |> 
  separate(X1,sep="\\|",into = c("c1","c2","c3","c4")) |> 
  select(-c("c1","c2","c4")) |> 
  rename(JGI_ID=c3,MMETSP_ID=X2,Percent_Identity=X3,Annotations=X10) |> 
  select(MMETSP_ID,JGI_ID,Percent_Identity,Annotations) |> 
  distinct() 

write_csv(blast,"proteomics/output_blasts/aurli_to_aurli/jgiblastannotations.csv")