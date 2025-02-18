library(pacman)
pacman::p_load(CHNOSZ,tidyverse,Biostrings)

# This script takes as input the blast of the Aurli proteome to Stramenopiles.
# and outputs annotations (Aurliprot_conserved_Stram.csv) based on
# whether it maps to Hondaea, Stramenopiles, Labys, or only to Aurli


#BLAST RESULTS
blast <- read_table("proteomics/output_blasts/all_aurli_to_stramenopiles/Aurliprot_Stram.out",
                    col_names = FALSE) |> 
  mutate(AN=paste(X10,X11,X12,X13,X14,X15)) |> 
  mutate(AN=gsub(" NA","",AN)) |> 
  select(-c(X10:X15,X9))

colnames(blast) <- c("qseqid","sseqid","pident","length","mismatch","gapopen",
                "evalue","bitscore","stitle")


blast <- blast |> 
  group_by(qseqid) |> 
  slice_min(order_by = evalue,n = 5) |> 
  slice_max(order_by = pident,n = 5) |> 
  mutate(Hondaea=ifelse(str_detect(stitle,"Hondaea"),TRUE,FALSE)) |> 
  mutate(Aurli=ifelse(str_detect(stitle,"Aurantio"),TRUE,FALSE)) 

Aurli <- blast|> 
  filter(Aurli==TRUE) |> 
  slice_min(order_by = evalue,n = 1) |> 
  slice_max(order_by = pident,n = 1)|> 
  slice_max(order_by = bitscore,n = 1)|> 
  slice_max(order_by = length,n = 1)|> 
  slice_min(order_by = mismatch,n = 1)|> 
  slice_min(order_by = gapopen,n = 1)|> 
  slice_sample(n=1)
HondF <- blast|> 
  filter(Aurli==FALSE) |> 
  filter(Hondaea==FALSE)|> 
  slice_min(order_by = evalue,n = 1) |> 
  slice_max(order_by = pident,n = 1)|> 
  slice_max(order_by = bitscore,n = 1)|> 
  slice_max(order_by = length,n = 1)|> 
  slice_min(order_by = mismatch,n = 1)|> 
  slice_min(order_by = gapopen,n = 1)|> 
  slice_sample(n=1)
HondT <- blast|> 
  filter(Hondaea==TRUE) |> 
  filter(Aurli==FALSE) |>
  slice_min(order_by = evalue,n = 1) |> 
  slice_max(order_by = pident,n = 1)|> 
  slice_max(order_by = bitscore,n = 1)|> 
  slice_max(order_by = length,n = 1)|> 
  slice_min(order_by = mismatch,n = 1)|> 
  slice_min(order_by = gapopen,n = 1)|> 
  slice_sample(n=1)


full_Table <- full_join(HondT,HondF,by="qseqid",suffix = c("_hondaea", "_stramenopiles")) |> 
  mutate(group=ifelse(is.na(Hondaea_hondaea)&!is.na(Hondaea_stramenopiles),"Stra_only",
                      ifelse(!is.na(Hondaea_hondaea)&is.na(Hondaea_stramenopiles),"Laby_only",
                             ifelse(!is.na(Hondaea_hondaea)&!is.na(Hondaea_stramenopiles),"both","error"))))

#Option 1: Aurli only
#No hits in the blast.
fastafile <- readAAStringSet("proteomics/input_fasta/Aurli1_aa.fa")
seq_name <- names(fastafile)

aurli_only <- seq_name[!seq_name %in% full_Table$qseqid]
#3190 in the proteome. 71 detected. 3 Significant

#Option 2:
#Protein found in Hondaea but not Stramenopiles. Laby only
# To do: Assume it is unique to Labys.
laby_only <- full_Table |> filter(group=="Laby_only") |> select(qseqid) |> pull()
#1290 in the proteome. 157 detected. 27 significant.

#Option 3:
#Protein not found in Hondaea but found in other Stramenopiles
# To do: Assume it was lost/not annotated in Hondaea, and it may be a common protein.
# Look also orthogrops in other labys.
not_hondaea <- full_Table |> filter(group=="Stra_only") |> select(qseqid) |> pull()
#2890 in the proteome. 910 detected. 3 Significant.


#Option 4:
#Protein found both in Hondaea and in other Stramenopiles
# To do: Substract %identity, E-value, bit-score (compare differences between Hondaea and Stramenopile
#with the highest similarity).

full_Table |> filter(group=="both") |> select(qseqid) |> pull()
#2445. Significant 457
both_table <- full_Table |> filter(group=="both")

# If difference is large (significant). Plot distribution of differences in bitscores.
#Assume protein is not well conserved in Labys maybe because it does a different function in Labys.
result <- both_table |> 
  mutate(diff_Tscore=bitscore_hondaea-bitscore_stramenopiles) |> 
  mutate(diff_E_Value=evalue_hondaea-evalue_stramenopiles) |> 
  mutate(diff_P_Ident=pident_hondaea-pident_stramenopiles) 

x <- quantile(result$diff_Tscore,probs=c(0.85,0.90,0.95,0.975))

result <- result |> 
  mutate(perc_diff_Tscore=ifelse(diff_Tscore<x[1],"0%-85%",
                                 ifelse(diff_Tscore<x[2]&diff_Tscore>=x[1],"85%-90%",
                                        ifelse(diff_Tscore<x[3]&diff_Tscore>=x[2],"90%-95%",
                                               ifelse(diff_Tscore<x[4]&diff_Tscore>=x[3],"95%-97.5%",
                                                      ifelse(diff_Tscore>=x[4],"97.5%-100%"))))))

ggplot(result, aes(x=diff_Tscore)) + 
  # geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
  #                binwidth=1,
  #                colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")+
  theme_minimal()+
  xlab("Aurli/Hondaea - Aurli/Stramen Total Score")+
  geom_vline(xintercept =0)+
  geom_vline(xintercept =x[1])+
  geom_vline(xintercept =x[2])+
  geom_vline(xintercept =x[3])+# Overlay with transparent density plot
  geom_vline(xintercept =x[4])# Overlay with transparent density plot



result <- result |> select(qseqid,diff_Tscore,diff_E_Value,diff_P_Ident,perc_diff_Tscore) |> 
  arrange(perc_diff_Tscore) |> 
  mutate(Group="Both")

a <- cbind(aurli_only,"aurli_only")
b <- cbind(laby_only,"laby_only")
c <- cbind(not_hondaea,"Stramenop_not_hondaea")
d <- rbind(a,b,c)
colnames(d) <- c("qseqid","Group")
d <- d |> 
  as.data.frame() |> 
  mutate(diff_Tscore=NA,
         diff_E_Value=NA,diff_P_Ident=NA,perc_diff_Tscore=NA) |> 
  select(qseqid,diff_Tscore,diff_E_Value,diff_P_Ident,perc_diff_Tscore,Group)

result1 <- rbind(d,result)

result1 |> filter(perc_diff_Tscore=="97.5%-100%") |> dim()

write.csv(result1,file = "proteomics/input_anno/Aurliprot_conserved_Stram.csv")
