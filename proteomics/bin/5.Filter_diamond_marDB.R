library(pacman)
pacman::p_load(CHNOSZ,tidyverse,Biostrings)

# This script takes as input the blast of the JGI and MMETSP proteome to Mariana's DB.
# and outputs annotations (Aurliprot_conserved_Stram.csv) based on
# whether it maps to Hondaea, Stramenopiles, Labys, or only to Aurli

setwd("/gpfs/projects/CollierGroup/agilgomez/projects/laby_proteomics/")

#BLAST RESULTS
 db <- "JGI_to_marDB.out"
#db <-  "mmetsp_to_marDB.out"
blast <- read_table(paste0("proteomics/output_blasts/diamond_to_marDB/",db),
                    col_names = FALSE)|>
                    select(-X9)

colnames(blast) <- c("qseqid","sseqid","pident","length","mismatch","gapopen",
                "evalue","bitscore","stitle")

tax_groups <- read.csv("proteomics/output_blasts/diamond_to_marDB/stramenopile_groups.csv")
str_groups <- tax_groups |> filter(Stramenopile=="TRUE")
not_str_groups <- tax_groups |> filter(is.na(Stramenopile))
stramenop <- paste0(str_groups$Taxonomic_group,collapse="|")
not_stramenop <- paste0(not_str_groups$Taxonomic_group,collapse="|")

laby_groups <- read.csv("proteomics/output_blasts/diamond_to_marDB/laby_groups.csv")
labys <- paste0(laby_groups$Taxonomic_group,collapse="|")

blast <- blast |>
  mutate(Aurantio=ifelse(str_detect(stitle,"Aurantio"),TRUE,FALSE))|>
  mutate(laby=ifelse(str_detect(stitle,labys),TRUE,FALSE)) |>
  mutate(strameno=ifelse(str_detect(stitle,stramenop),TRUE,FALSE))|>
  mutate(non_strameno=ifelse(str_detect(stitle,not_stramenop),TRUE,FALSE))
 
blast_sum <- blast |> 
  group_by(qseqid) |>
  summarise(aurantio_s=sum(Aurantio),
            laby_s=sum(laby),
            strameno_s=sum(strameno),
            non_strameno_s=sum(non_strameno))|>
  mutate(category=ifelse(aurantio_s==laby_s & laby_s == strameno_s & non_strameno_s==0,"Aurantio_only",
        ifelse(laby_s == strameno_s & non_strameno_s==0,"Laby_only",
        ifelse(laby_s != strameno_s,"Laby_and_Strameno",
        ifelse(non_strameno_s==0,"Strameno_only",
        ifelse(non_strameno_s>0 & strameno_s>0,"Strameno_and_Nonstrameno","ERROR"))))))

blast_sum |>
select(qseqid,category) |> distinct() |> select(category) |> table()
blast1 <- left_join(blast,blast_sum,by="qseqid")

# Group 0: All things in the blast have one copy in Aurli (the query).

# Group 1: Aurantio only (only found in the Aurantio sp genomes). Aurantio=Laby=Stramenopiles
aurantio_only <- blast1 |> filter(category=="Aurantio_only") |> select(qseqid) |> pull() |> unique()
aurantio_only |> length() # 3623 JGI; 174 MMETSP

# Group 2: Laby only. Protein found in Labys but not in Stramenopiles. Laby=Stramenopiles.  
laby_only <- blast1 |> filter(category=="Laby_only") |> select(qseqid) |> pull() |> unique()
laby_only |> length() # 8672 JGI; 2846 MMETSP

#Group 3: Protein found in Labys and in other Stramenopiles. Laby!=Stramenopiles
# To do: Substract %identity, E-value, bit-score (compare differences between Laby and Stramenopile
#with the highest similarity).
strameno_laby <- blast1 |> filter(category=="Laby_and_Strameno") |> select(qseqid) |> pull() |> unique()
strameno_laby |> length() # 1359 JGI; 112 MMETSP

LabyT <- blast1|> 
  filter(category=="Laby_and_Strameno") |> 
  filter(str_detect(stitle,labys))|>
  group_by(qseqid)|>
  slice_min(order_by = evalue,n = 1) |> 
  slice_max(order_by = pident,n = 1)|> 
  slice_max(order_by = bitscore,n = 1)|> 
  slice_max(order_by = length,n = 1)|> 
  slice_min(order_by = mismatch,n = 1)|> 
  slice_min(order_by = gapopen,n = 1)|> 
  slice_sample(n=1)

LabyF <- blast1|> 
  filter(category=="Laby_and_Strameno") |> 
  filter(!str_detect(stitle,labys))|>
  group_by(qseqid)|>
  slice_min(order_by = evalue,n = 1) |> 
  slice_max(order_by = pident,n = 1)|> 
  slice_max(order_by = bitscore,n = 1)|> 
  slice_max(order_by = length,n = 1)|> 
  slice_min(order_by = mismatch,n = 1)|> 
  slice_min(order_by = gapopen,n = 1)|> 
  slice_sample(n=1)

Laby_str_Table <- full_join(LabyT,LabyF,
by="qseqid",suffix = c("_laby", "_stramenopiles"))

result1_strameno_laby <- Laby_str_Table |> 
  mutate(diff_Tscore=bitscore_laby-bitscore_stramenopiles) |> 
  mutate(diff_E_Value=evalue_laby-evalue_stramenopiles) |> 
  mutate(diff_P_Ident=pident_laby-pident_stramenopiles) 

x <- quantile(result1_strameno_laby$diff_Tscore,probs=c(0.85,0.90,0.95,0.975))

result1_strameno_laby <- result1_strameno_laby |> 
  mutate(perc_diff_Tscore=ifelse(diff_Tscore<x[1],"0%-85%",
                                 ifelse(diff_Tscore<x[2]&diff_Tscore>=x[1],"85%-90%",
                                        ifelse(diff_Tscore<x[3]&diff_Tscore>=x[2],"90%-95%",
                                               ifelse(diff_Tscore<x[4]&diff_Tscore>=x[3],"95%-97.5%",
                                                      ifelse(diff_Tscore>=x[4],"97.5%-100%"))))))

result_l_s <- result1_strameno_laby |> select(qseqid,diff_Tscore,diff_E_Value,diff_P_Ident,perc_diff_Tscore) |> 
  arrange(perc_diff_Tscore) 

ggplot(result1_strameno_laby, aes(x=diff_Tscore)) + 
  # geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
  #                binwidth=1,
  #                colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")+
  theme_minimal()+
  xlab("Aurli/Labys - Aurli/Stramen Total Score")+
  geom_vline(xintercept =0)+
  geom_vline(xintercept =x[1])+
  geom_vline(xintercept =x[2])+
  geom_vline(xintercept =x[3])+# Overlay with transparent density plot
  geom_vline(xintercept =x[4])# Overlay with transparent density plot

#Group 4: Protein found in Stramenopiles_only but not in Laby
strameno_only <- blast1 |> filter(category=="Strameno_only") |> select(qseqid) |> pull()|> unique()
strameno_only |> length() # 0

# Group 5: Protein found in Stramenopiles and other Non-stramenopiles
# To do: Substract %identity, E-value, bit-score (compare differences between Stramenopile and non-Stramenopile
#with the highest similarity).
non_strameno <- blast1 |> filter(category=="Strameno_and_Nonstrameno") |> select(qseqid) |> pull() |> unique()
non_strameno |> length() # 1139 JGI; 158 MMETSP

strT <- blast1|> 
  filter(category=="Strameno_and_Nonstrameno") |> 
  filter(str_detect(stitle,labys))|>
  group_by(qseqid)|>
  slice_min(order_by = evalue,n = 1) |> 
  slice_max(order_by = pident,n = 1)|> 
  slice_max(order_by = bitscore,n = 1)|> 
  slice_max(order_by = length,n = 1)|> 
  slice_min(order_by = mismatch,n = 1)|> 
  slice_min(order_by = gapopen,n = 1)|> 
  slice_sample(n=1)

strF <- blast1|> 
  filter(category=="Strameno_and_Nonstrameno") |> 
  filter(str_detect(stitle,not_stramenop))|>
  group_by(qseqid)|>
  slice_min(order_by = evalue,n = 1) |> 
  slice_max(order_by = pident,n = 1)|> 
  slice_max(order_by = bitscore,n = 1)|> 
  slice_max(order_by = length,n = 1)|> 
  slice_min(order_by = mismatch,n = 1)|> 
  slice_min(order_by = gapopen,n = 1)|> 
  slice_sample(n=1)


Str_nStr_Table <- full_join(strT,strF,
by="qseqid",suffix = c("_stramenopiles", "_nonstramenopiles"))

result1_strameno_nonstra <- Str_nStr_Table |> 
  mutate(diff_Tscore=bitscore_stramenopiles-bitscore_nonstramenopiles) |> 
  mutate(diff_E_Value=evalue_stramenopiles-evalue_nonstramenopiles) |> 
  mutate(diff_P_Ident=pident_stramenopiles-pident_nonstramenopiles) 

x <- quantile(result1_strameno_nonstra$diff_Tscore,probs=c(0.85,0.90,0.95,0.975))

result1_strameno_nonstra <- result1_strameno_nonstra |> 
  mutate(perc_diff_Tscore=ifelse(diff_Tscore<x[1],"0%-85%",
                                 ifelse(diff_Tscore<x[2]&diff_Tscore>=x[1],"85%-90%",
                                        ifelse(diff_Tscore<x[3]&diff_Tscore>=x[2],"90%-95%",
                                               ifelse(diff_Tscore<x[4]&diff_Tscore>=x[3],"95%-97.5%",
                                                      ifelse(diff_Tscore>=x[4],"97.5%-100%"))))))

# If difference is large (significant). Plot distribution of differences in bitscores.
#Assume protein is not well conserved in Labys maybe because it does a different function in Labys.

ggplot(result1_strameno_nonstra, aes(x=diff_Tscore)) + 
  # geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
  #                binwidth=1,
  #                colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")+
  theme_minimal()+
  xlab("Aurli/Laby - Aurli/Non-Strameno Total Score")+
  geom_vline(xintercept =0)+
  geom_vline(xintercept =x[1])+
  geom_vline(xintercept =x[2])+
  geom_vline(xintercept =x[3])+# Overlay with transparent density plot
  geom_vline(xintercept =x[4])# Overlay with transparent density plot

result_s_ns <- result1_strameno_nonstra |> select(qseqid,diff_Tscore,diff_E_Value,diff_P_Ident,perc_diff_Tscore) |> 
  arrange(perc_diff_Tscore) 

####
# Combine results:
group_a <- cbind(aurantio_only,"aurantio_only")
group_b <- cbind(laby_only,"laby_only")
df <- rbind(group_a,group_b)
colnames(df) <- c("qseqid","Group")
df <- df |> 
  as.data.frame() |> 
  mutate(diff_Tscore=NA,
         diff_E_Value=NA,diff_P_Ident=NA,perc_diff_Tscore=NA) |> 
  select(qseqid,diff_Tscore,diff_E_Value,diff_P_Ident,perc_diff_Tscore,Group)

result_l_s <- result_l_s |>
mutate(Group="Laby_and_Strameno")

result_s_ns <- result_s_ns |>
mutate(Group="Laby_and_Non_Strameno")

result_all <- rbind(df,result_l_s,result_s_ns)

result_all |> filter(perc_diff_Tscore=="97.5%-100%") |> dim()

write.csv(result_all,file = paste0("proteomics/input_anno/Aurliprot_conserved_",db,"_Stram.csv"))
