library(pacman)
pacman::p_load(tidyverse,
               Biostrings)

# This script takes as input the blast of the JGI and MMETSP proteome to Mariana's DB.
# and outputs annotations (Aurliprot_conserved_Stram.csv) based on
# whether it maps to Hondaea, Stramenopiles, Labys, or only to Aurli

setwd("/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/")

#BLAST RESULTS
# database_name <- "JGI"
# seq_path <- 'proteomics/input_fasta/Aurli1_GeneCatalog_proteins_20120618.aa.fasta'
 database_name <- "mmetsp"
seq_path <- 'proteomics/input_fasta/mmetsp_uniprot.fa'
# database_name <-"mmetspENA"
#seq_path <- 'proteomics/input_fasta/ox87102.2020_06.faa'


fastafile <- readAAStringSet(seq_path)
seq_name <- names(fastafile)
sequence <- paste(fastafile)
df <- data.frame(seq_name, sequence) |>
  as_tibble() 
df$seq_name |> unique() |> length()

diam_colnames <- c("qseqid","sseqid","pident","length","mismatch","gapopen",
                "evalue","bitscore","stitle")
aurantio_only <- read_table(paste0("proteomics/output_blasts/diamond_to_marDB/",
                    database_name,"_to_Aurantio_only.out"),
                    col_names = FALSE)
colnames(aurantio_only) <- diam_colnames
laby_not_aurantio <- read_table(paste0("proteomics/output_blasts/diamond_to_marDB/",
                    database_name,"_to_Laby_notAurantio.out"),
                    col_names = FALSE)
colnames(laby_not_aurantio) <- diam_colnames
str_not_laby <- read_table(paste0("proteomics/output_blasts/diamond_to_marDB/",
                    database_name,"_to_str_not_laby.out"),
                    col_names = FALSE)
colnames(str_not_laby) <- diam_colnames
nonstr <- read_table(paste0("proteomics/output_blasts/diamond_to_marDB/",
                    database_name,"_to_nonstr_only.out"),
                    col_names = FALSE)
colnames(nonstr) <- diam_colnames

Aurli_sum <- aurantio_only |> 
  mutate(Aurli=ifelse(str_detect(stitle,"Aurli"),TRUE,FALSE))|>
  group_by(qseqid) |>
  summarise(total_Aurli=sum(Aurli))

Aurantio_not_Aurli_sum <- aurantio_only |> 
  mutate(not_Aurli=ifelse(str_detect(stitle,"Aurli"),FALSE,TRUE))|>
  group_by(qseqid) |>
  summarise(total_Aurantio_not_Aurli=sum(not_Aurli))

# Number of hits to Aurli
aurantio_only$qseqid |> unique() |> length()
Aurli_sum |> dim()

sum(Aurli_sum$total_Aurli ==50)
max(Aurli_sum$total_Aurli)
# 0 proteins only have hits to Aurli (JGI). 50 hits. Max 35 in Aurli
# 0  proteins only have hits to Aurli (mmetsp). 50 hits. Max 30 in Aurli
# 0 proteins only have hits to Aurli (mmetspENA). 50 hits. Max 33 in Aurli
Aurantio <- left_join(Aurli_sum,Aurantio_not_Aurli_sum,by="qseqid")
laby_not_aurantio$qseqid |> unique() |> length()

aurantio_only$qseqid |> unique() |> length() # 14816 JGI

# Labys not Aurantio:
Laby_sum <- laby_not_aurantio |> 
  group_by(qseqid) |>
  summarise(total_Laby_not_Aurantio=n())
laby_not_aurantio$qseqid |> unique() |> length() # 12088 JGI

Str_not_Laby_sum <- str_not_laby |> 
  group_by(qseqid) |>
  summarise(total_Str_not_laby=n())
str_not_laby$qseqid |> unique() |> length() # 9528 JGI

Nonstr_sum <- nonstr |> 
  group_by(qseqid) |>
  summarise(total_nonStr=n())
nonstr$qseqid |> unique() |> length() # 9821 JGI

Counts <- Aurantio |>
  left_join(Laby_sum,by="qseqid")  |>
  left_join(Str_not_Laby_sum,by="qseqid") |> 
  left_join(Nonstr_sum,by="qseqid")|>
  mutate(total_Laby_not_Aurantio=ifelse(is.na(total_Laby_not_Aurantio),0,total_Laby_not_Aurantio))|>
  mutate(total_Str_not_laby=ifelse(is.na(total_Str_not_laby),0,total_Str_not_laby))|>
  mutate(total_nonStr=ifelse(is.na(total_nonStr),0,total_nonStr))

blast_sum <- Counts |> 
  group_by(qseqid) |>
  mutate(Aurli_only= ifelse(total_Aurli>0 & total_Aurantio_not_Aurli==0 & total_Laby_not_Aurantio ==0 &
        total_Str_not_laby==0 & total_nonStr ==0,TRUE,FALSE))|>
  mutate(Aurantio_only= ifelse((total_Aurli > 0 | total_Aurantio_not_Aurli > 0) & total_Laby_not_Aurantio==0 & 
        total_Str_not_laby ==0 & total_nonStr==0,TRUE,FALSE))|>
  mutate(Laby_only= ifelse((total_Aurli > 0 | total_Aurantio_not_Aurli > 0 | total_Laby_not_Aurantio > 0) &
         total_Str_not_laby ==0 & total_nonStr==0,TRUE,FALSE))|>
  mutate(Str_only= ifelse((total_Aurli > 0 | total_Aurantio_not_Aurli > 0 | total_Laby_not_Aurantio > 0 |  total_Str_not_laby >0) &
         total_nonStr==0,TRUE,FALSE))|>
  mutate(category=ifelse(Aurli_only,"Aurli_only",
  ifelse(Aurantio_only,"Aurantio_only",
  ifelse(Laby_only,"Laby_only",
  ifelse(Str_only,"Str_only",
  ifelse(Aurli_only& !Aurantio_only & !Laby_only & Str_only,"Aurli_Stramenopile",
  ifelse(Aurli_only& !Aurantio_only & !Laby_only & !Str_only,"Aurli_NonStramenopile",
  ifelse(Aurli_only& Aurantio_only & !Laby_only & Str_only,"Aurantio_Stramenopile",
  ifelse(Aurli_only& Aurantio_only & !Laby_only & !Str_only,"Aurantio_NonStramenopile",
  ifelse(Aurli_only& Aurantio_only & Laby_only & Str_only,"Laby_Stramenopile",
  ifelse(Aurli_only& Aurantio_only & Laby_only & !Str_only,"Laby_NonStramenopile","Eukaryote")))))))))))
  
gr <- blast_sum$category |> table()
gr_p <- round(100*gr/sum(gr),3)
gr_x <- cbind(gr,gr_p)
write.table(gr_x,file = paste0("proteomics/input_anno/Aurliprot_conserved_",database_name,"_splitstats.txt"))

aurantio_only$db <- "Aurantio_only"
laby_not_aurantio$db <- "Laby_notAurantio"
str_not_laby$db <- "Str_notLaby"
nonstr$db <- "NonStr"
all_diamond <- rbind(aurantio_only,laby_not_aurantio,str_not_laby,nonstr)

best_diamond <- all_diamond |>
group_by(db,qseqid)|>
  slice_min(order_by = evalue,n = 1) |> 
  slice_max(order_by = pident,n = 1)|> 
  slice_max(order_by = bitscore,n = 1)|> 
  slice_max(order_by = length,n = 1)|> 
  slice_min(order_by = mismatch,n = 1)|> 
  slice_min(order_by = gapopen,n = 1)|> 
  slice_sample(n=1)

comb_table <- left_join(best_diamond,blast_sum,by="qseqid")

# Group 0: All things in the blast have one copy in Aurli (the query).

# Group 1: Aurantio only (only found in the Aurantio sp genomes).
aurli_only <- comb_table |> filter(category=="Aurli_only")
aurli_only$qseqid |> unique() |> length() # 1908 JGI; x MMETSP; x MMETSP ENA

aurantio_only <- comb_table |> filter(category=="Aurantio_only")
aurantio_only$qseqid |> unique() |> length()  # 768 JGI; x MMETSP; x MMETSP ENA

# Group 2: Laby only. Protein found in Labys but not in Stramenopiles.
laby_only <- comb_table |> filter(category=="Laby_only")
laby_only$qseqid |> unique() |> length()  # 2069 JGI; x MMETSP; x MMETSP ENA

# Group 3: Str_only
str_only <- comb_table |> filter(category=="Str_only")
str_only$qseqid |> unique() |> length() # 250 JGI; x MMETSP; x MMETSP ENA

StrnotLaby_T <- str_only |> 
  filter(db=="Str_notLaby")|>
  group_by(qseqid)|>
  slice_min(order_by = evalue,n = 1) |> 
  slice_max(order_by = pident,n = 1)|> 
  slice_max(order_by = bitscore,n = 1)|> 
  slice_max(order_by = length,n = 1)|> 
  slice_min(order_by = mismatch,n = 1)|> 
  slice_min(order_by = gapopen,n = 1)|> 
  slice_sample(n=1)

StrnotLaby_F <- str_only |> 
  filter(db!="Str_notLaby")|>
  group_by(qseqid)|>
  slice_min(order_by = evalue,n = 1) |> 
  slice_max(order_by = pident,n = 1)|> 
  slice_max(order_by = bitscore,n = 1)|> 
  slice_max(order_by = length,n = 1)|> 
  slice_min(order_by = mismatch,n = 1)|> 
  slice_min(order_by = gapopen,n = 1)|> 
  slice_sample(n=1)

dim(StrnotLaby_T)
dim(StrnotLaby_F)

StrnotLaby_table <- full_join(StrnotLaby_F,StrnotLaby_T,
by="qseqid",suffix = c("_StrnotLaby_F","_StrnotLaby_T"))

result1_strameno_TF <- StrnotLaby_table |> 
filter(!is.na(bitscore_StrnotLaby_F))|>
filter(!is.na(bitscore_StrnotLaby_T))|>
  mutate(diff_Tscore=bitscore_StrnotLaby_F-bitscore_StrnotLaby_T) |> 
  mutate(diff_E_Value=evalue_StrnotLaby_F-evalue_StrnotLaby_T) |> 
  mutate(diff_P_Ident=pident_StrnotLaby_F-pident_StrnotLaby_T) 

x <- quantile(result1_strameno_TF$diff_Tscore,probs=c(0.85,0.90,0.95,0.975))

result1_strameno_TF <- result1_strameno_TF |> 
  mutate(perc_diff_Tscore=ifelse(diff_Tscore<x[1],"0%-85%",
                                 ifelse(diff_Tscore<x[2]&diff_Tscore>=x[1],"85%-90%",
                                        ifelse(diff_Tscore<x[3]&diff_Tscore>=x[2],"90%-95%",
                                               ifelse(diff_Tscore<x[4]&diff_Tscore>=x[3],"95%-97.5%",
                                                      ifelse(diff_Tscore>=x[4],"97.5%-100%"))))))

result_l_s <- result1_strameno_TF |> select(qseqid,diff_Tscore,diff_E_Value,diff_P_Ident,perc_diff_Tscore) |> 
  arrange(perc_diff_Tscore) 

ggplot(result1_strameno_TF, aes(x=diff_Tscore)) + 
  # geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
  #                binwidth=1,
  #                colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")+
  theme_minimal()+
  xlab("Laby - Stramenopiles Total Score")+
  geom_vline(xintercept =0)+
  geom_vline(xintercept =x[1])+
  geom_vline(xintercept =x[2])+
  geom_vline(xintercept =x[3])+# Overlay with transparent density plot
  geom_vline(xintercept =x[4])# Overlay with transparent density plot


# Group 4: Eukaryotes
euk <- comb_table |> filter(category=="Eukaryote")
euk$qseqid |> unique() |> length() # 9821 JGI; x MMETSP; x MMETSP ENA

euk_T <- euk |> 
  filter(db=="NonStr"| db=="Str_notLaby")|>
  group_by(qseqid)|>
  slice_min(order_by = evalue,n = 1) |> 
  slice_max(order_by = pident,n = 1)|> 
  slice_max(order_by = bitscore,n = 1)|> 
  slice_max(order_by = length,n = 1)|> 
  slice_min(order_by = mismatch,n = 1)|> 
  slice_min(order_by = gapopen,n = 1)|> 
  slice_sample(n=1)

euk_F <- euk |> 
  filter(db!="NonStr" | db!="Str_notLaby")|>
  group_by(qseqid)|>
  slice_min(order_by = evalue,n = 1) |> 
  slice_max(order_by = pident,n = 1)|> 
  slice_max(order_by = bitscore,n = 1)|> 
  slice_max(order_by = length,n = 1)|> 
  slice_min(order_by = mismatch,n = 1)|> 
  slice_min(order_by = gapopen,n = 1)|> 
  slice_sample(n=1)

dim(euk_T)
dim(euk_F)

euk_table <- full_join(euk_F,euk_T,
by="qseqid",suffix = c("_euk_F","_euk_T"))

result1_euk_TF <- euk_table |> 
filter(!is.na(bitscore_euk_F))|>
filter(!is.na(bitscore_euk_T))|>
  mutate(diff_Tscore=bitscore_euk_F-bitscore_euk_T) |> 
  mutate(diff_E_Value=evalue_euk_F-evalue_euk_T) |> 
  mutate(diff_P_Ident=pident_euk_F-pident_euk_T) 

x <- quantile(result1_euk_TF$diff_Tscore,probs=c(0.85,0.90,0.95,0.975))

result1_euk_TF <- result1_euk_TF |> 
  mutate(perc_diff_Tscore=ifelse(diff_Tscore<x[1],"0%-85%",
                                 ifelse(diff_Tscore<x[2]&diff_Tscore>=x[1],"85%-90%",
                                        ifelse(diff_Tscore<x[3]&diff_Tscore>=x[2],"90%-95%",
                                               ifelse(diff_Tscore<x[4]&diff_Tscore>=x[3],"95%-97.5%",
                                                      ifelse(diff_Tscore>=x[4],"97.5%-100%"))))))

result_l_ns <- result1_euk_TF |> select(qseqid,diff_Tscore,diff_E_Value,diff_P_Ident,perc_diff_Tscore) |> 
  arrange(perc_diff_Tscore) 

ggplot(result1_euk_TF, aes(x=diff_Tscore)) + 
  # geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
  #                binwidth=1,
  #                colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")+
  theme_minimal()+
  xlab("Laby - Eukaryote Total Score")+
  geom_vline(xintercept =0)+
  geom_vline(xintercept =x[1])+
  geom_vline(xintercept =x[2])+
  geom_vline(xintercept =x[3])+# Overlay with transparent density plot
  geom_vline(xintercept =x[4])# Overlay with transparent density plot


####
# Combine results:
group_a <- cbind(unique(aurli_only$qseqid),"aurli_only")
group_b <- cbind(unique(aurantio_only$qseqid),"aurantio_only")
group_c <- cbind(unique(laby_only$qseqid),"laby_only")
df <- rbind(group_a,group_b,group_c)
colnames(df) <- c("qseqid","Group")
df <- df |> 
  as.data.frame() |> 
  mutate(diff_Tscore=NA,
         diff_E_Value=NA,diff_P_Ident=NA,perc_diff_Tscore=NA) |> 
  select(qseqid,diff_Tscore,diff_E_Value,diff_P_Ident,perc_diff_Tscore,Group)

result_l_s$qseqid |> unique() |> length()
result_l_ns$qseqid |> unique() |> length()

result_l_s <- result_l_s |>
mutate(Group="Stramenopile")

result_l_ns <- result_l_ns |>
mutate(Group="Eukaryote")

result_all <- rbind(df,result_l_s,result_l_ns)

write.csv(result_all,file = paste0("proteomics/input_anno/Aurliprot_conserved_",database_name,"_split.csv"),row.names=FALSE)
