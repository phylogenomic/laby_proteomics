library(rentrez)
library(tidyverse)

# Replace empty lines with >
#awk '!NF{$0=">"}1' file

line <- read_lines("proteomics/results/inter_intra_blast/short_files_ncbiweb/alltables1.txt") 

line <- ifelse(str_starts(line,"Alignments:"),">",line)
line <- ifelse(str_starts(line,"Sequences producing"),">",line)
line <- ifelse(str_starts(line,"Alignments:"),">",line)
line <- ifelse(str_starts(line,"Sequence ID:"),">",line)
line <- ifelse(str_starts(line,"Range"),">",line)
line <- ifelse(str_starts(line,"Score:"),">",line)
line <- ifelse(str_starts(line,"Method:"),">",line)
line <- ifelse(str_starts(line,">"),">",line)
line <- ifelse(str_detect(line,"Scientific"),">",line)
line <- ifelse(str_detect(line,"Name"),">",line)
line <- ifelse(str_detect(line,"No significant"),"NS;NS;NS;NS;NS;NS;NS;NS;NS",line)
line <- str_squish(line)
line <- line[line!="--"&line!=">"]


columns_names <- c("Description","Taxid","Max_Score","Total_Score","Query_cover","E_Value",
"Percent_Identity","Accesion_Length","Accession")

vec <- c()
for(i in 1:length(line)){
  if(str_starts(line[i],"Query #")){
    c <- str_split(line[i]," ",n=8)[[1]]
    vec[i] <- paste0(c[3],";",c[8])}
  else if(str_starts(line[i],"NS;NS;NS;NS")){
    vec[i] <- line[i]
  }
  else{
    c <- str_split(line[i]," ",n=30)[[1]]
    L <- length(c)
    to <- L-8
    name1 <- paste0(c[1:to],collapse = " ")
    to2 <- to+1
    rest <- c[to2:L]
    rest1 <- paste0(rest,collapse = ";")
    vec[i] <- paste0(name1,";",rest1)
  }
}

write_lines(vec,"alltables2.txt")


df <- data.frame(matrix(ncol = 2, nrow = 0))
for (i in 1:length(vec)){
  if (str_starts(vec[i],"jgi")){
    df[i,1] <- vec[i]
    jgi <- vec[i]
  } else{
    df[i,1] <- jgi
    df[i,2] <- vec[i]
  }
}

df1 <- df |> 
  filter(!is.na(X2)) |> 
  separate(X1,into=c("Query","Length"),";") |> 
  separate(X2,into=columns_names,";") 




#Option 1:
#Protein not found in Stramenopiles (other than Aurli).
# To do: Assume it is unique to Aurli.
G1 <- df1 |> 
  filter(Description=="NS")
aurli_only <- G1$Query
#3

# Filter 10 per unique query.
unique <- df1 |> 
  filter(Description!="NS") |> 
  mutate(Max_Score=as.numeric(Max_Score),
         Total_Score=as.numeric(Total_Score),
         Accesion_Length=as.numeric(Accesion_Length),
         E_Value=as.numeric(E_Value),
         Percent_Identity=as.numeric(Percent_Identity)) |> 
  group_by(Query) |> 
  slice_min(order_by = E_Value,n = 5) |> 
  slice_max(order_by = Percent_Identity,n = 5) |> 
  mutate(Hondaea=ifelse(Taxid=="2315210","TRUE","FALSE")) 


HondF <- unique|> 
  filter(Taxid!="2315210")|> 
  slice_min(order_by = E_Value,n = 1) |> 
  slice_max(order_by = Percent_Identity,n = 1)|> 
  slice_max(order_by = Total_Score,n = 1)|> 
  slice_max(order_by = Query_cover,n = 1)|> 
  slice_max(order_by = Max_Score,n = 1)|> 
  slice_max(order_by = Accesion_Length,n = 1)|> 
  slice_sample(n=1)
HondT <- unique|> 
  filter(Taxid=="2315210")|> 
  slice_min(order_by = E_Value,n = 1) |> 
  slice_max(order_by = Percent_Identity,n = 1)|> 
  slice_max(order_by = Total_Score,n = 1)|> 
  slice_max(order_by = Query_cover,n = 1)|> 
  slice_max(order_by = Max_Score,n = 1)|> 
  slice_max(order_by = Accesion_Length,n = 1) |> 
  slice_sample(n=1)

full_Table <- full_join(HondT,HondF,by="Query",suffix = c("_hondaea", "_stramenopiles")) |>  
  mutate(group=ifelse(is.na(Taxid_hondaea)&!is.na(Taxid_stramenopiles),"Stra_only",
                      ifelse(!is.na(Taxid_hondaea)&is.na(Taxid_stramenopiles),"Laby_only",
                             ifelse(!is.na(Taxid_hondaea)&!is.na(Taxid_stramenopiles),"both","error"))))


#Option 2:
#Protein found in Hondaea but not Stramenopiles.
# To do: Assume it is unique to Labys.
laby_only <- full_Table |> filter(group=="Laby_only") |> select(Query) |> pull()
#27

#Option 3:
#Protein not found in Hondaea but found in other Stramenopiles
# To do: Assume it was lost/not annotated in Hondaea, and it may be a common protein.
# Look also orthogrops in other labys.
not_hondaea <- full_Table |> filter(group=="Stra_only") |> select(Query) |> pull()
#3

#Option 4:
#Protein found both in Hondaea and in other Stramenopiles
# To do: Substract %identity, E-value, bit-score (compare differences between Hondaea and Stramenopile
#with the highest similarity).

full_Table |> filter(group=="both") |> select(Query) |> pull()
#457
both_table <- full_Table |> filter(group=="both")

# If difference is large (significant). Plot distribution of differences in bitscores.
  #Assume protein is not well conserved in Labys maybe because it does a different function in Labys.
result <- both_table |> 
  mutate(diff_Tscore=Total_Score_hondaea-Total_Score_stramenopiles) |> 
  mutate(diff_E_Value=E_Value_hondaea-E_Value_stramenopiles) |> 
  mutate(diff_P_Ident=Percent_Identity_hondaea-Percent_Identity_stramenopiles) 

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



result <- result |> select(Query,diff_Tscore,diff_E_Value,diff_P_Ident,perc_diff_Tscore) |> 
  arrange(perc_diff_Tscore) |> 
  mutate(Group="Both")

a <- cbind(aurli_only,"aurli_only")
b <- cbind(laby_only,"laby_only")
c <- cbind(not_hondaea,"Stramenop_not_hondaea")
d <- rbind(a,b,c)
colnames(d) <- c("Query","Group")
d <- d |> 
  as.data.frame() |> 
  mutate(diff_Tscore=NA,
         diff_E_Value=NA,diff_P_Ident=NA,perc_diff_Tscore=NA) |> 
  select(Query,diff_Tscore,diff_E_Value,diff_P_Ident,perc_diff_Tscore,Group)

result1 <- rbind(d,result)

write.csv(result1,file = "proteomics/results/inter_intra_blast/percen_table.csv")
