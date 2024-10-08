####################################
# Differential enrichment analysis for Proteomics data.
# author: "Alejandro Gil-Gomez"
library(pacman)
pacman::p_load(tidyverse, Biostrings, FactoMineR,
               data.table, plotly, DEP, sva,
               SummarizedExperiment, ComplexHeatmap,
               patchwork, factoextra,
               ggrepel, RColorBrewer,
               ggVennDiagram, grid, here,ggplotify,
               ggvenn)

contr <- "t0" 
#contr <- "tx" 

fold_threshold <- 1.0001 #1.2, 
inverse_fold <- 1 / fold_threshold
p_value <- 0.1 # 0.05

col_gradient <- c("#ffffd4","#fed98e","#fe9929","#d95f0e","#993404")

file_set <- c("mmetsp","jgi","all_det","all_sig")

## Load data

# JGI SET
filename <- "rawnormalized_revised.csv"

raw_data <- read_csv(file.path("proteomics/input_data", filename))
colnames(raw_data) <- c("name_ID", "Fasta_headers",
                        "RAW.T0-1", "RAW.T0-2", "RAW.T0-3",
                        "RAW.T2-1", "RAW.T2-2", "RAW.T2-3",
                        "RAW.T4-1", "RAW.T4-2", "RAW.T4-3",
                        "RAW.T6-1", "RAW.T6-2", "RAW.T6-3",
                        "RAW.T8-1", "RAW.T8-2", "RAW.T8-3")

data <- raw_data

experimental_design <- data.frame(label = c("T0-1", "T0-2", "T0-3",
                                            "T2-1", "T2-2", "T2-3",
                                            "T4-1", "T4-2", "T4-3",
                                            "T6-1", "T6-2", "T6-3",
                                            "T8-1", "T8-2", "T8-3"),
                                  condition = c("T0", "T0", "T0",
                                                "T2", "T2", "T2",
                                                "T4", "T4", "T4",
                                                "T6", "T6", "T6",
                                                "T8", "T8", "T8"),
                                  replicate = rep(1:3, 5))

data_jgi <- data |>
  separate(name_ID, sep = "\\|", into = letters[1:4]) |>
  select(-c("a", "b", "d")) |>
  rename("name_ID" = "c") |>
  select(-Fasta_headers)

data_jgi$set <- "jgi"

# MMETSP SET
filename <- "rawnormalized.csv"

raw_data <- read_csv(file.path("proteomics/input_data", filename))
colnames(raw_data) <- c("name_ID", "Fasta_headers",
                        "RAW.T0-1", "RAW.T0-2", "RAW.T0-3",
                        "RAW.T2-1", "RAW.T2-2", "RAW.T2-3",
                        "RAW.T4-1", "RAW.T4-2", "RAW.T4-3",
                        "RAW.T6-1", "RAW.T6-2", "RAW.T6-3",
                        "RAW.T8-1", "RAW.T8-2", "RAW.T8-3")

data <- raw_data

experimental_design <- data.frame(label = c("T0-1", "T0-2", "T0-3",
                                            "T2-1", "T2-2", "T2-3",
                                            "T4-1", "T4-2", "T4-3",
                                            "T6-1", "T6-2", "T6-3",
                                            "T8-1", "T8-2", "T8-3"),
                                  condition = c("T0", "T0", "T0",
                                                "T2", "T2", "T2",
                                                "T4", "T4", "T4",
                                                "T6", "T6", "T6",
                                                "T8", "T8", "T8"),
                                  replicate = rep(1:3, 5))

# MMETSPs IDs.
data_mmetsp <- data |> 
  select(-Fasta_headers) 
data_mmetsp$set <- "mmetsp"

# Combine all sets.
raw_sets <- list(data_mmetsp = data_mmetsp,
                 data_jgi = data_jgi)

raw_sets <-  raw_sets |>
  map(ungroup) |>
  map(mutate, Protein_id = paste0(name_ID, "-", set)) |> 
  map(mutate, name=name_ID) |> 
  map(mutate, ID=name_ID) |> 
  map(select,Protein_id,name_ID,name,ID,starts_with("RAW"),set)

raw_sets |> 
  map(dim) # Total number of proteins


# Are there any duplicate entries?
raw_sets |>
  map(ungroup) |>
  map(group_by, name_ID) |>
  map(summarize, frequency = n()) |>
  map(arrange, desc(frequency)) |>
  map(filter, frequency > 1)

# Rawsets done!

# DEP pipeline:
DEP_sets <- raw_sets |>
  map(make_unique, "name", "ID", delim = ";")

# How many proteins were detected?
DEP_sets |>
  map(dim) 
# 3298 mmetsp, 3006 jgi

DEP_sets |>
  map(colnames) 

# Remove batch effect using SVA
DEP_sets_numbers <- DEP_sets |>
  map(select, starts_with("RAW"))
batch_corrected <- DEP_sets_numbers |>
  map(ComBat, as.character(rep(1:3, 5)))
for (i in seq_along(batch_corrected)){
  colnames(batch_corrected[[i]]) <- gsub("RAW","BATCH_CORRECT",colnames(batch_corrected[[i]]))
}

DEP_sets <- DEP_sets |>
  map(select, c(1:4))
DEP_batch_corrected <- map2(DEP_sets, batch_corrected, cbind)
DEP_batch_corrected |>
  map(dim)

# Save RAW_BATCH CORRECTED as csv
raw_batch <- map2(raw_sets,DEP_batch_corrected,left_join,
                  by="name",
                  suffix = c("", "_n")) |> 
  map(select,!ends_with("_n"))

for (i in seq_along(raw_batch)) {
  raw_batch[[i]] |>
    write_csv(file.path("proteomics/DEP_results_IndependentSets", paste0("raw_batch",file_set[i],".",contr,".csv")))
}

#Columns with the replicates/time points must be in the same order than SE.
batch_c_names <- DEP_batch_corrected |> 
  map(select,starts_with("BATCH")) |> 
  map(colnames)

seq <- 1:length(DEP_batch_corrected[[1]])

lfq_columns <- seq[colnames(DEP_batch_corrected[[1]]) %in% batch_c_names[[1]]]

DEP_bc_se <- DEP_batch_corrected |>
  map(make_se, lfq_columns, experimental_design)
# Data diagnostics.
DEP_bc_se |>
  map(plot_frequency)
data_filt <- DEP_bc_se |>
  map(filter_missval, thr = 0)
data_filt |> map(plot_numbers)
data_filt |> map(plot_coverage)
# Normalize datasets
data_norm <- data_filt |>
  map(normalize_vsn)
map2(data_filt, data_norm, plot_normalization)
data_norm |> map(plot_normalization)

# Detect and remove missing values (if any)
# data_filt |> map(plot_missval)
# data_filt |> map(plot_detect)
# data_imp <- data_norm |> map(impute, fun = "MinProb", q = 0.01)
# map2(data_norm, data_imp, plot_imputation)

# Test differences with respect to T0
# data_diff <- data_norm |>
#  map(test_diff, type = "control", control = "T0")
# data_diff_all_contrasts <- data_norm |>
#  map(test_diff, type = "all")
# dep <- data_diff |> map(add_rejections, alpha = p_value, lfc = log2(fold_threshold))
# dep_get_results <- dep |> map(get_results)
# dep_get_wide <- dep |> map(get_df_wide)

#Test contrasts:
if (contr == "t0") {
  data_diff_manual <- data_norm |>
    map(test_diff, type = "manual", test = c("T2_vs_T0","T4_vs_T0","T6_vs_T0","T8_vs_T0"
    ))
} else if (contr == "tx") {
  data_diff_manual <- data_norm |>
    map(test_diff, type = "manual", test = c("T2_vs_T0","T4_vs_T0","T6_vs_T0","T8_vs_T0",
                                             "T4_vs_T2",
                                             "T6_vs_T4",
                                             "T8_vs_T6"
    ))
} else {
  print("error with contrasts")
}

dep <- data_diff_manual |>
  map(add_rejections,alpha = p_value, lfc = log2(fold_threshold))

#save(dep,file="proteomics/dep_sets.R")

dep_get_results <- dep |> 
  map(get_results)
row_names <- dep_get_results|>
  map(rownames) |>
  map(as.data.frame) |>
  map(rename, rownumber = ".x[[i]]")
dep_get_results <- map2(dep_get_results, row_names, cbind)

# Add mapping: name_anno
mapping1 <- read_csv("proteomics/input_anno/jgimmetsp_result_may24.csv") |> 
  rename(name=MMETSP_ID) |> 
  rename(name_anno=JGI_ID)

dep_get_results_an <- list()
for (i in 1:length(dep_get_results)){
  if(names(dep_get_results)[[i]]=="data_mmetsp"){
    # map mmetsp to jgi
    
    dep_get_results_an[[i]] <- left_join(dep_get_results[[i]],mapping1,by="name") |>
      mutate(name_anno=ifelse(is.na(name_anno),name,name_anno)) |> 
      select(name,name_anno,rownumber,3:24)
    
  }else if(names(dep_get_results)[i]=="data_jgi"){
    # do nothing
    dep_get_results_an[[i]] <- dep_get_results[[i]]|> 
      mutate(name_anno=as.character(name)) |> 
      select(name,name_anno,rownumber,3:24)
  }
}

# Venn diagrams:
# Venn by size in nVennR. Venn normal is ggVennDiagram.
# Venn of detected. GGVenndiagram
detected_name <- dep_get_results_an |>
  map(select, name, name_anno) 

mmetsp_det <- detected_name[[1]] |>
  select(name_anno) |> arrange(desc(name_anno)) |> 
  pull()
  
mmetsp_det_M <- mmetsp_det[str_detect(mmetsp_det,"A")]
mmetsp_det_MAP <- mmetsp_det[str_detect(mmetsp_det,"A",negate=TRUE)]
jgi_det <- detected_name[[2]] |> select(name_anno) |> pull()
setdiff(mmetsp_det,jgi_det)

y <- list("MMETSP only" = mmetsp_det_M,
          "MMETSP mapped" = mmetsp_det_MAP,
          "JGI" = jgi_det)
# 4D Venn diagram
venn_detected <- ggVennDiagram(y, lwd = 0.8, lty = 1) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  ggtitle("Detected proteins")

# Venn of significants
significant_name <- dep_get_results_an |>
  map(filter, significant == TRUE) |>
  map(select, name,name_anno) 

mmetsp_sig <- significant_name[[1]]|>
  select(name_anno) |> arrange(desc(name_anno)) |> 
  pull()
# MMETSP proteins that didn't map to JGI
mmetsp_sig_M <- mmetsp_sig[str_detect(mmetsp_sig,"A")]
# MMETSP proteins that map to JGI
mmetsp_sig_MAP <- mmetsp_sig[str_detect(mmetsp_sig,"A",negate=TRUE)] 
jgi_sig <- significant_name[[2]]|> select(name_anno) |> pull()

setdiff(mmetsp_sig_MAP,jgi_sig)

y <- list("MMETSP only" = mmetsp_sig_M,
          "MMETSP mapped" = mmetsp_sig_MAP,
          "JGI" = jgi_sig)
# 4D Venn diagram
venn_significant <- ggVennDiagram(y, lwd = 0.8, lty = 1) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  ggtitle("Significant proteins")

# Both
y <- list("MMETSP detected" = mmetsp_det,
          "JGI detected" = jgi_det,
          "MMETSP significant" = mmetsp_sig,
          "JGI significant" = jgi_sig)

y <- list("MMETSP only detected" = mmetsp_det_M,
          "MMETSP mapped detected" = mmetsp_det_MAP,
          "JGI detected" = jgi_det,
          "MMETSP only sig" = mmetsp_sig_M,
          "MMETSP mapped sig" = mmetsp_sig_MAP,
          "JGI sig" = jgi_sig)
# 4D Venn diagram
venn_both <- ggVennDiagram(y, lwd = 0.8, lty = 1) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  ggtitle("Proteins")

venn_plot <- wrap_plots(venn_detected, venn_significant,venn_both) +
  plot_annotation(tag_levels = "A", tag_suffix = ".") &
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 15, hjust = 0, vjust = 0))

ggsave(
  filename = paste0("proteomics/img_IndependentSets/indset.fig2_t0_vennplot.pdf"),
  plot = venn_plot,
  width = 15,
  height = 13
)

ggsave(
  filename = paste0("proteomics/img_IndependentSets/indset.fig2_t0_vennplot.png"),
  plot = venn_plot,
  width = 15,
  height = 13
)


# Create all_det set:
# For proteins significant in both take MMETSP version,
# If there is a dispute (significant/not significant M vs J), take the sign one.
venn_detected

mmetsp_det_names <- detected_name[[1]] |> 
  filter(name_anno %in% mmetsp_det) # Detected in MMETSP
mmetsp_sig_names <- significant_name[[1]] |> 
  filter(name_anno %in% mmetsp_sig) # Significant in MMETSP
mmetsp_notsig_names <- setdiff(mmetsp_det_names$name_anno,
                               mmetsp_sig_names$name_anno) # Not significant in MMETSP
jgi_det <- detected_name[[2]] |> 
  select(name_anno) |> pull() # Detected in JGI
jgi_sig <- significant_name[[2]] |> 
  select(name_anno) |> pull() # Significant in JGI
jgi_notsig <- setdiff(jgi_det,jgi_sig) # Not significant in JGI

# Possibilities for all detected set (+ means significant)
# M J Final_detected
# + + Detected_both Get_M 
msig_jsig <- intersect(mmetsp_sig_names$name_anno,jgi_sig) 
msig_jsig |> length() # 301 sig M and sig J
msig_jsig_names <- detected_name[[1]] |> 
  filter(name_anno %in% msig_jsig) |> 
  select(name) |> pull()
msig_jsig_names |> length() # 307 MMETSP names match name_anno
#msig_jsig_names <- msig_jsig_names[msig_jsig_names %in% mmetsp_sig_names$name]
#msig_jsig_names |> length() # 302 MMETSP names
dep1_msig_jsig_names <- dep[[1]][rownames(dep[[1]])%in%msig_jsig_names,]

# + - Detected_both Get_M 
msig_jnotsig <- intersect(mmetsp_sig_names$name_anno,jgi_notsig) 
msig_jnotsig |> length() # 150 sig M and not J
msig_jnotsig_names <- detected_name[[1]] |> 
  filter(name_anno %in% msig_jnotsig) |> 
  select(name) |> pull()
#msig_jnotsig_names <- msig_jnotsig_names[msig_jnotsig_names %in% mmetsp_sig_names$name]
msig_jnotsig_names |> length() # 152
dep1_msig_jnotsig_names <- dep[[1]][rownames(dep[[1]])%in%msig_jnotsig_names,]

# - + Detected_both Get_J.
mnotsig_jsig <- intersect(mmetsp_notsig_names,jgi_sig) 
mnotsig_jsig |> length() # 52 not M and sig J
mnotsig_jsig_names <- mnotsig_jsig 
mnotsig_jsig_names |> length() #52
dep2_mnotsig_jsig_names <- dep[[2]][rownames(dep[[2]])%in%mnotsig_jsig_names,]

# - - Detected_both Get_M. 2040
mnotsig_jnotsig <- intersect(mmetsp_notsig_names,jgi_notsig) # 2014 not M and not J
mnotsig_jnotsig |> length() # 2014 not M and not J
mnotsig_jnotsig_names <- detected_name[[1]] |> 
  filter(name_anno %in% mnotsig_jnotsig) |> 
  select(name) |> pull() 
mnotsig_jnotsig_names |> length() #2040
dep1_mnotsig_jnotsig_names <- dep[[1]][rownames(dep[[1]])%in%mnotsig_jnotsig_names,]

# Detect None. Absent:
mmetsp_detected_notjgi <- setdiff(mmetsp_det_names$name_anno,jgi_det)
mmetsp_detected_notjgi |> length() # 737 Detected M not in J
jgi_detected_notmmetsp <- setdiff(jgi_det,mmetsp_det_names$name_anno)
jgi_detected_notmmetsp |> length()  # 489 Detected J not in M

# + 0 Detected_M Get_M.
msig_jnotdet <- intersect(mmetsp_sig_names$name_anno,
                          mmetsp_detected_notjgi) 
msig_jnotdet|> length()# 77 Sig M and not found in J
msig_jnotdet_names <- detected_name[[1]] |> 
  filter(name_anno %in% msig_jnotdet) |> 
  select(name) |> pull()
#msig_jnotdet_names <- msig_jnotdet_names[msig_jnotdet_names %in% mmetsp_sig_names$name]
msig_jnotdet_names |> length() #81
dep1_msig_jnotdet_names <- dep[[1]][rownames(dep[[1]])%in%msig_jnotdet_names,]

# 0 + Detected_J Get_J.
mnotdet_jsig <- intersect(jgi_detected_notmmetsp,jgi_sig) 
mnotdet_jsig |> length() # 41 sig J and not found in M
mnotdet_jsig_names <- detected_name[[2]] |> 
  filter(name_anno %in% mnotdet_jsig) |> 
  select(name) |> pull()
mnotdet_jsig_names |> length() #41
dep2_mnotdet_jsig_names <- dep[[2]][rownames(dep[[2]])%in%mnotdet_jsig_names,]

# - 0 Detected_M Get_M.
mnotsig_jnotdet <- intersect(mmetsp_notsig_names,
                             mmetsp_detected_notjgi)
mnotsig_jnotdet |> length() # 660 notS in M and not found in J
mnotsig_jnotdet_names <- detected_name[[1]] |> 
  filter(name_anno %in% mnotsig_jnotdet) |> 
  select(name) |> pull()
mnotsig_jnotdet_names |> length() #662
dep1_mnotsig_jnotdet_names <- dep[[1]][rownames(dep[[1]])%in%mnotsig_jnotdet_names,]

# 0 - Detected_J Get_J.
mnotdet_jnotsig <- intersect(jgi_detected_notmmetsp,
                             jgi_notsig) 
mnotdet_jnotsig |> length() # 448 notS in J and not found in M
mnotdet_jnotsig_names <- detected_name[[2]] |> 
  filter(name_anno %in% mnotdet_jnotsig) |> 
  select(name) |> pull()#448
mnotdet_jnotsig_names |> length()
dep2_mnotdet_jnotsig_names <- dep[[2]][rownames(dep[[2]])%in%mnotdet_jnotsig_names,]

# 0 0 Not found. 0
intersect(jgi_detected_notmmetsp,
          mmetsp_detected_notjgi) # 0 not found

# Make all_det following the rules:
# For proteins significant in both take MMETSP version,
# If there is a dispute (significant/not M/J),
# take the one that is significant.
c(msig_jsig_names, # 302 sig M and sig J
  msig_jnotsig_names, # 151 sig M and notS J
  mnotsig_jsig_names,# 52 notS M and sig J
  mnotsig_jnotsig_names, # 2040 notS M and notS J
  msig_jnotdet_names, # 77 Sig M and not found in J
  mnotdet_jsig_names, # 41 sig J and not found in M
  mnotsig_jnotdet_names, # 662 notS in M and not found in J
  mnotdet_jnotsig_names)|>  # 448 notS in J and not found in M
  unique() |> 
  length() # 3783

# Make all_det following the rules:
# For proteins significant in both take MMETSP version,
# If there is a dispute (significant/not M/J), take the one that is significant.
dep_alldet <- rbind(dep1_msig_jsig_names, # SET M: 307 significant both
                    dep1_msig_jnotsig_names, # SET M: 152 significant M
                    dep2_mnotsig_jsig_names, # SET J: 52 significant in J
                    dep1_mnotsig_jnotsig_names, # SET M: 2040 not significant
                    dep1_msig_jnotdet_names, # SET M: 81 significant M
                    dep2_mnotdet_jsig_names, # SET J: 41 significant J
                    dep1_mnotsig_jnotdet_names, # SET M: 662 not significant
                    dep2_mnotdet_jnotsig_names) # SET J: 448 not significant
# 3783

dep[[3]] <- dep_alldet
names(dep)[3] <- "data_alldet"

# Create all_sig set:
venn_significant
mmetsp_sig_names <- significant_name[[1]] |> 
  filter(name_anno %in% mmetsp_sig) |> 
  select(name) |> pull() # 530
jgisig_only_ids <- setdiff(jgi_sig,mmetsp_sig_MAP) # 93

dep1_sig <- dep[[1]][rownames(dep[[1]])%in%mmetsp_sig_names,]
dep2_sig <- dep[[2]][rownames(dep[[2]])%in%jgisig_only_ids,]
dep[[4]] <- rbind(dep1_sig,dep2_sig)
names(dep)[4] <- "data_allsig"

dep[[4]] |> dim() # 623
significant_name

# Count total significant:
dep |>
  map(get_results) |>
  map(filter, significant == TRUE) |>
  map(summarise, n())
# 530, 394, 623, 623

# 
x <- dep |>
map(get_results) |>
  map(filter, significant == TRUE) 

sig_in_alldet <- x[[3]]$name
expected_sig <- c(rownames(dep1_msig_jsig_names), # SET M: 307 significant both
                  rownames(dep1_msig_jnotsig_names), # SET M: 152 significant M
                  rownames(dep2_mnotsig_jsig_names), # SET J: 52 significant in M
                  rownames(dep1_msig_jnotdet_names), # SET M: 81 significant M
                  rownames(dep2_mnotdet_jsig_names) # SET J: 41 significant J
                  ) # 633
sig_miss <- setdiff(expected_sig,sig_in_alldet)
setdiff(rownames(dep1_msig_jnotdet_names),sig_in_alldet)

# 5 first in 1st set. "A0A6S8BS59" "A0A6S8C0L7" "A0A6S8C3V9" "A0A6S8E961" "A0A6S8F057"  
# 1 in 2nd set. "A0A6S8G6P9"
# 4 in last set. "A0A6S8B3X8" "A0A6S8DAW0" "A0A6S8DMD4" "A0A6S8F133"
  
dep_get_wide <- dep |> map(get_df_wide)

# DEP
dep_get_results <- dep |> 
  map(get_results)
row_names <- dep_get_results|>
  map(rownames) |>
  map(as.data.frame) |>
  map(rename, rownumber = ".x[[i]]")
dep_get_results <- map2(dep_get_results, row_names, cbind)
dep_get_results |> 
  map(colnames) |> 
  map(length)

#Merge results with annotations:
# Combine result columns (raw, batch, get_results and get_wide)
if (contr == "t0") {
  results_base <- map2(raw_sets,DEP_batch_corrected,left_join,by="name",
                       suffix = c("", "_n")) |> 
    map(select,!ends_with("_n"))
  
  results1 <- dep_get_results |> 
    map2(dep_get_wide,left_join,by="name",
         suffix = c("", "_n")) |> 
    map(select,!ends_with("_n")) |> 
    map(select,Protein_id,name_ID,name,rownumber,
        #starts_with("RAW."),starts_with("BATCH_"),
        T0_1,T0_2,T0_3,
        T2_1,T2_2,T2_3,
        T4_1,T4_2,T4_3,
        T6_1,T6_2,T6_3,
        T8_1,T8_2,T8_3,
        ends_with("_ratio"),
        ends_with("_centered"),
        ends_with("_p.val"),
        ends_with("_p.adj"),
        ends_with("significant"),
        starts_with("T2_vs_T0_"),
        starts_with("T4_vs_T0_"),
        starts_with("T6_vs_T0_"),
        starts_with("T8_vs_T0_")) |> 
    map(~ mutate(.x, across(starts_with("T") & ends_with("ratio"), ~ 2^., .names = "FC_{.col}"))) |> 
    map(~ .x  |> rename_with(~ str_remove(., "_ratio"), starts_with("FC") & ends_with("_ratio")))
  
} else if (contr == "tx") {
  results_base <- map2(raw_sets,DEP_batch_corrected,left_join,by="name",
                       suffix = c("", "_n")) |> 
    map(select,!ends_with("_n"))
  
  results1 <- dep_get_results |> 
    map2(dep_get_wide,left_join,by="name",
         suffix = c("", "_n")) |> 
    map(select,!ends_with("_n")) |> 
    map(select,Protein_id,name_ID,name,rownumber,
        starts_with("RAW."),starts_with("BATCH_"),
        T0_1,T0_2,T0_3,
        T2_1,T2_2,T2_3,
        T4_1,T4_2,T4_3,
        T6_1,T6_2,T6_3,
        T8_1,T8_2,T8_3,
        ends_with("_ratio"),
        ends_with("_centered"),
        ends_with("_p.val"),
        ends_with("_p.adj"),
        ends_with("significant"),
        starts_with("T2_vs_T0_"),
        starts_with("T4_vs_T0_"),
        starts_with("T6_vs_T0_"),
        starts_with("T8_vs_T0_")) |> 
    map(~ mutate(.x, across(starts_with("T") & ends_with("ratio"), ~ 2^., .names = "FC_{.col}"))) |> 
    map(~ .x  |> rename_with(~ str_remove(., "_ratio"), starts_with("FC") & ends_with("_ratio")))
  
} else {
  # Code to execute if 'set' is neither "t0" nor "tx"
  # Add your code here
}

results1 |> 
  map(colnames) |> 
  map(length)

# Columns named the following:
# RAW.T#_#: Raw intensities. Not present in results1
# BATCH_CORRECT.T#_#: Intensities after batch correction (Combat::SVA)
# T#_#: (log2 of batch corrected columns after DEP::VSN normalization)
# T#_vs_T0_ratio: (ratio of log2(intensities)). Ex: T2_vs_T0_ratio=mean(T2_1,T2_2,T2_3)-mean(T0_1,T0_2,T0_3). It is centered in 0.
# T#_centered: Average log2 FC scaled protein wise. Ex: T0_centered=mean(T0_1,T0_2,T0_3)-mean(T0_1,T0_2,T0_3,T2_1,T2_2,T2_3,T4_1,T4_2,T4_3,T6_1,T6_2,T6_3,T8_1,T8_2,T8_3)
# T#_vs_T0_p.val
# T#_vs_T0_p.adj
# T#_vs_T0_significant
# T#_vs_T0_diff: diff is another way to calculate the ratio. # ratio and diff are linearly correlated but they are not the same.
# T#_vs_T0_CI.L: CI of diff.
# T#_vs_T0_CI.R: CI of diff.


# DEP results:
if (contr == "t0") {
  results2 <- results1 |> map(mutate,
                              T2v0.reg = ifelse(FC_T2_vs_T0 > fold_threshold & T2_vs_T0_p.adj <= p_value, "UP",
                                                ifelse(FC_T2_vs_T0 < inverse_fold & T2_vs_T0_p.adj <= p_value, "DOWN", "NO")),
                              T4v0.reg = ifelse(FC_T4_vs_T0 > fold_threshold & T4_vs_T0_p.adj <= p_value, "UP",
                                                ifelse(FC_T4_vs_T0 < inverse_fold & T4_vs_T0_p.adj <= p_value, "DOWN", "NO")),
                              T6v0.reg = ifelse(FC_T6_vs_T0 > fold_threshold & T6_vs_T0_p.adj <= p_value, "UP",
                                                ifelse(FC_T6_vs_T0 < inverse_fold & T6_vs_T0_p.adj <= p_value, "DOWN", "NO")),
                              T8v0.reg = ifelse(FC_T8_vs_T0 > fold_threshold & T8_vs_T0_p.adj <= p_value, "UP",
                                                ifelse(FC_T8_vs_T0 < inverse_fold & T8_vs_T0_p.adj <= p_value, "DOWN", "NO")))
  # Expression greater than the median.
  median_t0 <- results2 |>
    map(select, Protein_id,T0_1,T0_2,T0_3) |>
    map(pivot_longer, names_to = "rep", values_to = "values", 2:4) |>
    map(group_by, Protein_id) |>
    map(summarise, mean_t0 = mean(values)) |>
    map(select, mean_t0) |>
    map(pull) |>
    map(median)
  
  for (i in seq_along(results2)){
    results2[[i]] <- results2[[i]] |> 
      mutate(above_median_exp_t0 = T0_1 > median_t0[[i]] | T0_2 > median_t0[[i]] | T0_3 > median_t0[[i]])
  }
  
  results2 <- results2 |> 
    map(mutate,
        UP_AND468 = ifelse(T4v0.reg == "UP" & T6v0.reg == "UP" & T8v0.reg == "UP",
                           TRUE, FALSE)) |>
    map(mutate,
        UP_OR468 = ifelse(T4v0.reg == "UP" | T6v0.reg == "UP" | T8v0.reg == "UP",
                          TRUE, FALSE)) |>
    map(mutate,
        DOWN_AND468 = ifelse(T4v0.reg == "DOWN" & T6v0.reg == "DOWN" & T8v0.reg == "DOWN",
                             TRUE, FALSE)) |>
    map(mutate,
        DOWN_OR468 = ifelse(T4v0.reg == "DOWN" | T6v0.reg == "DOWN" | T8v0.reg == "DOWN",
                            TRUE, FALSE))  |> 
    map(mutate,
        Increasing = ifelse(T8_vs_T0_diff > T6_vs_T0_diff &
                              T6_vs_T0_diff > T4_vs_T0_diff &
                              T4_vs_T0_diff > T2_vs_T0_diff,
                            TRUE, FALSE)) |>
    map(ungroup)
  
} else if (contr == "tx") {
  results2 <- results1 |> map(mutate,
                              T2v0.reg = ifelse(FC_T2_vs_T0 > fold_threshold & T2_vs_T0_p.adj <= p_value, "UP",
                                                ifelse(FC_T2_vs_T0 < inverse_fold & T2_vs_T0_p.adj <= p_value, "DOWN", "NO")),
                              T4v0.reg = ifelse(FC_T4_vs_T0 > fold_threshold & T4_vs_T0_p.adj <= p_value, "UP",
                                                ifelse(FC_T4_vs_T0 < inverse_fold & T4_vs_T0_p.adj <= p_value, "DOWN", "NO")),
                              T6v0.reg = ifelse(FC_T6_vs_T0 > fold_threshold & T6_vs_T0_p.adj <= p_value, "UP",
                                                ifelse(FC_T6_vs_T0 < inverse_fold & T6_vs_T0_p.adj <= p_value, "DOWN", "NO")),
                              T8v0.reg = ifelse(FC_T8_vs_T0 > fold_threshold & T8_vs_T0_p.adj <= p_value, "UP",
                                                ifelse(FC_T8_vs_T0 < inverse_fold & T8_vs_T0_p.adj <= p_value, "DOWN", "NO")),
                              T4v2.reg = ifelse(FC_T4_vs_T2 > fold_threshold & T4_vs_T2_p.adj <= p_value, "UP",
                                                ifelse(FC_T4_vs_T2 < inverse_fold & T4_vs_T2_p.adj <= p_value, "DOWN", "NO")),
                              T6v4.reg = ifelse(FC_T6_vs_T4 > fold_threshold & T6_vs_T4_p.adj <= p_value, "UP",
                                                ifelse(FC_T6_vs_T4 < inverse_fold & T6_vs_T4_p.adj <= p_value, "DOWN", "NO")),
                              T8v6.reg = ifelse(FC_T8_vs_T6 > fold_threshold & T8_vs_T6_p.adj <= p_value, "UP",
                                                ifelse(FC_T8_vs_T6 < inverse_fold & T8_vs_T6_p.adj <= p_value, "DOWN", "NO")))
  
  # Expression greater than the median.
  median_t0 <- results2 |>
    map(select, Protein_id,T0_1,T0_2,T0_3) |>
    map(pivot_longer, names_to = "rep", values_to = "values", 2:4) |>
    map(group_by, Protein_id) |>
    map(summarise, mean_t0 = mean(values)) |>
    map(select, mean_t0) |>
    map(pull) |>
    map(median)
  
  for (i in seq_along(results2)){
    results2[[i]] <- results2[[i]] |> 
      mutate(above_median_exp_t0 = T0_1 > median_t0[[i]] | T0_2 > median_t0[[i]] | T0_3 > median_t0[[i]])
  }
  
  results2 <- results2 |> 
    map(mutate,
        UP_AND468 = ifelse(T4v0.reg == "UP" & T6v0.reg == "UP" & T8v0.reg == "UP",
                           TRUE, FALSE)) |>
    map(mutate,
        UP_OR468 = ifelse(T4v0.reg == "UP" | T6v0.reg == "UP" | T8v0.reg == "UP",
                          TRUE, FALSE)) |>
    map(mutate,
        DOWN_AND468 = ifelse(T4v0.reg == "DOWN" & T6v0.reg == "DOWN" & T8v0.reg == "DOWN",
                             TRUE, FALSE)) |>
    map(mutate,
        DOWN_OR468 = ifelse(T4v0.reg == "DOWN" | T6v0.reg == "DOWN" | T8v0.reg == "DOWN",
                            TRUE, FALSE)) |>
    map(mutate,
        UP_T6T4_OR_T4T2 = ifelse(T6v4.reg == "UP" | T4v2.reg == "UP",
                                 TRUE, FALSE)) |> 
    map(mutate,
        Increasing = ifelse(T8_vs_T0_diff > T6_vs_T0_diff &
                              T6_vs_T0_diff > T4_vs_T0_diff &
                              T4_vs_T0_diff > T2_vs_T0_diff,
                            TRUE, FALSE)) |>
    map(ungroup)
} else {
  # Code to execute if 'set' is neither "t0" nor "tx"
  # Add your code here
}

results2 |> map(colnames) |> map(length)

# Add mapping: name_anno
results3 <- list()
for (i in 1:length(results2)){
  if(names(results2)[[i]]=="data_mmetsp"){
    # map mmetsp to jgi
    results3[[i]] <- left_join(results2[[i]],mapping1,by="name") |>
      mutate(name_anno=ifelse(is.na(name_anno),name,name_anno)) |> 
      select(Protein_id,name_ID,name,name_anno,rownumber,5:67)
    
  }else if(names(results2)[i]=="data_jgi"){
    # do nothing
    results3[[i]] <- results2[[i]]|> 
      mutate(name_anno=as.character(name)) |> 
      select(Protein_id,name_ID,name,name_anno,rownumber,5:67)
  }else if(names(results2)[i]=="data_alldet"){
    df_mmetsp <- results2[[i]][str_detect(results2[[i]]$name,"A"),]
    df_jgi <- results2[[i]][str_detect(results2[[i]]$name,"A",negate=TRUE),]
    
    df_mmetsp1 <- left_join(df_mmetsp,mapping1,by="name") |>
      mutate(name_anno=ifelse(is.na(name_anno),name,name_anno)) |> 
      select(Protein_id,name_ID,name,name_anno,rownumber,5:67)
    
    df_jgi1 <- df_jgi|> 
      mutate(name_anno=as.character(name)) |> 
      select(Protein_id,name_ID,name,name_anno,rownumber,5:67)
    
    results3[[i]] <- rbind(df_mmetsp1,df_jgi1)
  }else if(names(results2)[i]=="data_allsig"){
    df_mmetsp <- results2[[i]][str_detect(results2[[i]]$name,"A"),]
    df_jgi <- results2[[i]][str_detect(results2[[i]]$name,"A",negate=TRUE),]
    
    df_mmetsp1 <- left_join(df_mmetsp,mapping1,by="name") |>
      mutate(name_anno=ifelse(is.na(name_anno),name,name_anno)) |> 
      select(Protein_id,name_ID,name,name_anno,rownumber,5:67)
    
    df_jgi1 <- df_jgi|> 
      mutate(name_anno=as.character(name)) |> 
      select(Protein_id,name_ID,name,name_anno,rownumber,5:67)
    
    results3[[i]] <- rbind(df_mmetsp1,df_jgi1)
  }
}
names(results3) <- c("data_mmetsp","data_jgi",
                     "data_alldet","data_allsig")
results3 |> map(colnames) |> map(length)


# Add RAW and BC
results4 <- list()
for (i in 1:length(results3)){
  if(names(results3)[[i]]=="data_mmetsp"){
    # map mmetsp to jgi
    results4[[i]] <- left_join(results3[[i]],raw_batch[[i]],
                               by="name",suffix = c("", "_n")) |> 
    select(!ends_with("_n"))
    
  }else if(names(results3)[i]=="data_jgi"){
    # do nothing
    results4[[i]] <- left_join(results3[[i]],raw_batch[[i]],
                               by="name",suffix = c("", "_n")) |> 
      select(!ends_with("_n"))
    
  }else if(names(results3)[i]=="data_alldet"){
    df_mmetsp <- results3[[i]][str_detect(results3[[i]]$name,"A"),]
    df_jgi <- results3[[i]][str_detect(results3[[i]]$name,"A",negate=TRUE),]
    
    df_mmetsp1 <- left_join(df_mmetsp,raw_batch[[1]],
                            by="name",suffix = c("", "_n")) |> 
      select(!ends_with("_n"))
    
    df_jgi1 <- left_join(df_jgi,raw_batch[[2]],
                         by="name",suffix = c("", "_n")) |> 
      select(!ends_with("_n"))
    
    results4[[i]] <- rbind(df_mmetsp1,df_jgi1)
  }else if(names(results3)[i]=="data_allsig"){
    df_mmetsp <- results3[[i]][str_detect(results3[[i]]$name,"A"),]
    df_jgi <- results3[[i]][str_detect(results3[[i]]$name,"A",negate=TRUE),]
    
    df_mmetsp1 <- left_join(df_mmetsp,raw_batch[[1]],
                            by="name",suffix = c("", "_n")) |> 
      select(!ends_with("_n"))
    
    df_jgi1 <- left_join(df_jgi,raw_batch[[2]],
                         by="name",suffix = c("", "_n")) |> 
      select(!ends_with("_n"))
    
    results4[[i]] <- rbind(df_mmetsp1,df_jgi1)
  }
}
results4 |> map(colnames) |> map(length)



############
#### MERGE ANNOTATIONS WITH EXTERNAL DATA.
ext_anno <- read_csv("proteomics/input_anno/all_anno_combined.csv") |> 
  mutate(name_anno=as.character(name_anno))

results_anno <- results4 |> 
  map(left_join,ext_anno,by="name_anno")

results_anno[[1]]$table <- "mmetsp"
results_anno[[2]]$table <- "jgi"
results_anno[[3]]$table <- "alldet"
results_anno[[4]]$table <- "allsig"

results_anno |> map(colnames)

# Annotations added!


#Venn Diagrams Significant
setnames <- c("MMETSP","JGI","ALLDET","ALLSIG")
which_set <- 3

if (contr == "t0") {
  t2v0_sig <- results_anno|>
    map(ungroup) |>
    map(filter, T2v0.reg == "DOWN" | T2v0.reg == "UP") |>
    map(select, name_anno) |>
    map(pull)
  t4v0_sig <- results_anno |>
    map(ungroup) |>
    map(filter, T4v0.reg == "DOWN" | T4v0.reg == "UP") |>
    map(select, name_anno) |>
    map(pull)
  t6v0_sig <- results_anno |>
    map(ungroup) |>
    map(filter, T6v0.reg == "DOWN" | T6v0.reg == "UP") |>
    map(select, name_anno) |>
    map(pull)
  t8v0_sig <- results_anno |>
    map(ungroup) |>
    map(filter, T8v0.reg == "DOWN" | T8v0.reg == "UP") |>
    map(select, name_anno) |>
    map(pull)
  
  tx_t0 <- list("T2 vs T0" = t2v0_sig[[which_set]],
                "T4 vs T0" = t4v0_sig[[which_set]],
                "T6 vs T0" = t6v0_sig[[which_set]],
                "T8 vs T0" = t8v0_sig[[which_set]])
  
  # 4D Venn diagram
  venn_sig <- ggVennDiagram(tx_t0, lwd = 0.8, lty = 1) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    ggtitle(paste(setnames[which_set]," significant proteins"))
  
  #Upregulated
  t2v0_up <- results_anno |>
    map(ungroup) |>
    map(filter, T2v0.reg == "UP") |>
    map(select, name_anno) |>
    map(pull)
  t4v0_up <- results_anno |>
    map(ungroup) |>
    map(filter, T4v0.reg == "UP") |>
    map(select, name_anno) |>
    map(pull)
  t6v0_up <- results_anno |>
    map(ungroup) |>
    map(filter, T6v0.reg == "UP") |>
    map(select, name_anno) |>
    map(pull)
  t8v0_up <- results_anno |>
    map(ungroup) |>
    map(filter, T8v0.reg == "UP") |>
    map(select, name_anno) |>
    map(pull)
  
  tx_t0 <- list("T2 vs T0" = t2v0_up[[which_set]],
                "T4 vs T0" = t4v0_up[[which_set]],
                "T6 vs T0" = t6v0_up[[which_set]],
                "T8 vs T0" = t8v0_up[[which_set]])
  
  # 4D Venn diagram
  venn_up <- ggVennDiagram(tx_t0, lwd = 0.8, lty = 1) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    ggtitle(paste(setnames[which_set]," upregulated proteins"))
  
  #Downregulated
  t2v0_down <- results_anno |>
    map(ungroup) |>
    map(filter, T2v0.reg == "DOWN") |>
    map(select, name_anno) |>
    map(pull)
  t4v0_down <- results_anno |>
    map(ungroup) |>
    map(filter, T4v0.reg == "DOWN") |>
    map(select, name_anno) |>
    map(pull)
  t6v0_down <- results_anno |>
    map(ungroup) |>
    map(filter, T6v0.reg == "DOWN") |>
    map(select, name_anno) |>
    map(pull)
  t8v0_down <- results_anno |>
    map(ungroup) |>
    map(filter, T8v0.reg == "DOWN") |>
    map(select, name_anno) |>
    map(pull)
  tx_t0 <- list("T2 vs T0" = t2v0_down[[which_set]],
                "T4 vs T0" = t4v0_down[[which_set]],
                "T6 vs T0" = t6v0_down[[which_set]],
                "T8 vs T0" = t8v0_down[[which_set]])
  
  # 4D Venn diagram
  venn_down <- ggVennDiagram(tx_t0, lwd = 0.8, lty = 1) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    ggtitle(paste(setnames[which_set]," downregulated proteins"))
  
  layout <- "
  ###AAAAA###
  ###AAAAA###
  BBBBB#CCCCC
  BBBBB#CCCCC
  "
  
  fig2.1 <- wrap_plots(venn_sig,
                       venn_up,venn_down)+
    plot_layout(design=layout)+
    plot_annotation(tag_levels = "A", tag_suffix = ".") &
    theme(plot.tag.position = c(0, 1),
          plot.tag = element_text(size = 15, hjust = 0, vjust = 0))
  
  ggsave(
    filename = paste0("proteomics/img_IndependentSets/indset.fig2_t0_",setnames[which_set],".pdf"),
    plot = fig2.1,
    width = 15,
    height = 13
  )
  
  ggsave(
    filename = paste0("proteomics/img_IndependentSets/indset.fig2_t0_",setnames[which_set],".png"),
    plot = fig2.1,
    width = 15,
    height = 13,dpi=500
  )
  
  
} else if (contr == "tx") {
  t2v0_sig <- results_anno|>
    map(ungroup) |>
    map(filter, T2v0.reg == "DOWN" | T2v0.reg == "UP") |>
    map(select, name_anno) |>
    map(pull)
  t4v0_sig <- results_anno |>
    map(ungroup) |>
    map(filter, T4v0.reg == "DOWN" | T4v0.reg == "UP") |>
    map(select, name_anno) |>
    map(pull)
  t6v0_sig <- results_anno |>
    map(ungroup) |>
    map(filter, T6v0.reg == "DOWN" | T6v0.reg == "UP") |>
    map(select, name_anno) |>
    map(pull)
  t8v0_sig <- results_anno |>
    map(ungroup) |>
    map(filter, T8v0.reg == "DOWN" | T8v0.reg == "UP") |>
    map(select, name_anno) |>
    map(pull)
  t4v2_sig <- results_anno |>
    map(ungroup) |>
    map(filter, T4v2.reg == "DOWN" | T4v2.reg == "UP") |>
    map(select, name_anno) |>
    map(pull)
  t6v4_sig <- results_anno |>
    map(ungroup) |>
    map(filter, T6v4.reg == "DOWN" | T6v4.reg == "UP") |>
    map(select, name_anno) |>
    map(pull)
  t8v6_sig <- results_anno |>
    map(ungroup) |>
    map(filter, T8v6.reg == "DOWN" | T8v6.reg == "UP") |>
    map(select, name_anno) |>
    map(pull)
  
  tx_t0 <- list("T2 vs T0" = t2v0_sig[[which_set]],
                "T4 vs T0" = t4v0_sig[[which_set]],
                "T6 vs T0" = t6v0_sig[[which_set]],
                "T8 vs T0" = t8v0_sig[[which_set]])
  
  tn2_tn1 <- list("T2 vs T0" = t2v0_sig[[which_set]],
                  "T4 vs T2" = t4v2_sig[[which_set]],
                  "T6 vs T4" = t6v4_sig[[which_set]],
                  "T8 vs T6" = t8v6_sig[[which_set]])
  
  # 4D Venn diagram
  venn_sig <- ggVennDiagram(tx_t0, lwd = 0.8, lty = 1) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    ggtitle(paste(setnames[which_set]," Merged significant proteins"))
  
  venn_sig2 <- ggVennDiagram(tn2_tn1, lwd = 0.8, lty = 1) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    ggtitle(paste(setnames[which_set]," Merged significant proteins"))
  
  #Upregulated
  t2v0_up <- results4 |>
    map(ungroup) |>
    map(filter, T2v0.reg == "UP") |>
    map(select, name_anno) |>
    map(pull)
  t4v0_up <- results4 |>
    map(ungroup) |>
    map(filter, T4v0.reg == "UP") |>
    map(select, name_anno) |>
    map(pull)
  t6v0_up <- results4 |>
    map(ungroup) |>
    map(filter, T6v0.reg == "UP") |>
    map(select, name_anno) |>
    map(pull)
  t8v0_up <- results4 |>
    map(ungroup) |>
    map(filter, T8v0.reg == "UP") |>
    map(select, name_anno) |>
    map(pull)
  t4v2_up <- results4 |>
    map(ungroup) |>
    map(filter, T4v2.reg == "UP") |>
    map(select, name_anno) |>
    map(pull)
  t6v4_up <- results4 |>
    map(ungroup) |>
    map(filter, T6v4.reg == "UP") |>
    map(select, name_anno) |>
    map(pull)
  t8v6_up <- results4 |>
    map(ungroup) |>
    map(filter, T8v6.reg == "UP") |>
    map(select, name_anno) |>
    map(pull)
  
  tx_t0 <- list("T2 vs T0" = t2v0_up[[which_set]],
                "T4 vs T0" = t4v0_up[[which_set]],
                "T6 vs T0" = t6v0_up[[which_set]],
                "T8 vs T0" = t8v0_up[[which_set]])
  
  tn2_tn1 <- list("T2 vs T0" = t2v0_up[[which_set]],
                  "T4 vs T2" = t4v2_up[[which_set]],
                  "T6 vs T4" = t6v4_up[[which_set]],
                  "T8 vs T6" = t8v6_up[[which_set]])
  
  
  # 4D Venn diagram
  venn_up <- ggVennDiagram(tx_t0, lwd = 0.8, lty = 1) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    ggtitle(paste(setnames[which_set]," upregulated proteins"))
  
  venn_up2 <- ggVennDiagram(tn2_tn1, lwd = 0.8, lty = 1) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    ggtitle(paste(setnames[which_set]," upregulated proteins"))
  
  #Downregulated
  t2v0_down <- results4 |>
    map(ungroup) |>
    map(filter, T2v0.reg == "DOWN") |>
    map(select, name_anno) |>
    map(pull)
  t4v0_down <- results4 |>
    map(ungroup) |>
    map(filter, T4v0.reg == "DOWN") |>
    map(select, name_anno) |>
    map(pull)
  t6v0_down <- results4 |>
    map(ungroup) |>
    map(filter, T6v0.reg == "DOWN") |>
    map(select, name_anno) |>
    map(pull)
  t8v0_down <- results4 |>
    map(ungroup) |>
    map(filter, T8v0.reg == "DOWN") |>
    map(select, name_anno) |>
    map(pull)
  t4v2_down <- results4 |>
    map(ungroup) |>
    map(filter, T4v2.reg == "DOWN") |>
    map(select, name_anno) |>
    map(pull)
  t6v4_down <- results4 |>
    map(ungroup) |>
    map(filter, T6v4.reg == "DOWN") |>
    map(select, name_anno) |>
    map(pull)
  t8v6_down <- results4 |>
    map(ungroup) |>
    map(filter, T8v6.reg == "DOWN") |>
    map(select, name_anno) |>
    map(pull)
  tx_t0 <- list("T2 vs T0" = t2v0_down[[which_set]],
                "T4 vs T0" = t4v0_down[[which_set]],
                "T6 vs T0" = t6v0_down[[which_set]],
                "T8 vs T0" = t8v0_down[[which_set]])
  
  tn2_tn1 <- list("T2 vs T0" = t2v0_down[[which_set]],
                  "T4 vs T2" = t4v2_down[[which_set]],
                  "T6 vs T4" = t6v4_down[[which_set]],
                  "T8 vs T6" = t8v6_down[[which_set]])
  
  # 4D Venn diagram
  venn_down <- ggVennDiagram(tx_t0, lwd = 0.8, lty = 1) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    ggtitle(paste(setnames[which_set]," downregulated proteins"))
  
  venn_down2 <- ggVennDiagram(tn2_tn1, lwd = 0.8, lty = 1) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    ggtitle("Merged downregulated proteins")
  
  layout <- "
  ###AAAAA###
  ###AAAAA###
  BBBBB#CCCCC
  BBBBB#CCCCC
  "
  
  fig2.2 <- wrap_plots(venn_sig2,
                       venn_up2,venn_down2)+
    plot_annotation(tag_levels = "A", tag_suffix = ".") &
    plot_layout(design=layout)+
    theme(plot.tag.position = c(0, 1),
          plot.tag = element_text(size = 15, hjust = 0, vjust = 0))
  
  ggsave(
    filename = paste0("proteomics/img_IndependentSets/indset.fig2_tx_",setnames[which_set],".pdf"),
    plot = fig2.2,
    width = 15,
    height = 13
  )
  
  ggsave(
    filename = paste0("proteomics/img_IndependentSets/indset.fig2_tx_",setnames[which_set],".png"),
    plot = fig2.2,
    width = 15,
    height = 13,dpi=500
  )
  
} else {
  # Code to execute if 'set' is neither "t0" nor "tx"
  # Add your code here
}

#################################
## DEP plots:

# PCA per replicate and condition
plot_pcas <- list()
title <- c("MMETSP set", "JGI set","AllDet set","AllSig set")
for (i in seq_along(dep)) {
  plot_pcas[[i]] <- plot_pca(dep[[i]], x = 1, y = 2,
                             n = dim(assay(dep[[i]]))[1], point_size = 4) +
    ggtitle(title[i]) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    theme(plot.title = element_text(size = 20))+
    scale_color_brewer(palette = "YlOrRd") +
    ylim(-9, 9)
}

# for Figure 4.
for (p in seq_along(plot_pcas)){
  fig4 <- plot_pcas[[p]] 
  ggsave(
    filename = paste0("proteomics/img_IndependentSets/indset.fig4.",file_set[p],".pdf"),
    plot = fig4,
    width = 7,
    height = 7
  )
  
  ggsave(
    filename = paste0("proteomics/img_IndependentSets/indset.fig4.",file_set[p],".png"),
    plot = fig4,
    width = 7,
    height = 7,dpi=500
  )
  
}


# Correlation plots
cor_plots <- dep |>
  map(plot_cor, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

assay1 <- dep |>
  map(SummarizedExperiment::assay) |>
  map(as.data.frame)

n <- assay1 |>
  map(rownames) |>
  map(as.data.frame) |>
  map(rename, rownumber = ".x[[i]]")
assay2 <- map2(assay1, n, cbind)

assay3 <- assay2 |>
  map(pivot_longer, names_to = "TimeRep",
      values_to = "Exp",
      1:15) |>
  map(separate, TimeRep, sep = "_", into = c("Time", "Rep")) |>
  map(mutate, TimeRep = paste0(Time, " ", Rep))

distribution_plots <- list()
for (i in seq_along(assay3)) {
  distribution_plots[[i]] <- assay3[[i]] |>
    ggplot(aes(x = Exp, y = TimeRep, fill = Time)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = "white") +
    theme_minimal() + scale_y_discrete(limits = rev) +
    theme(legend.position = "none") +
    xlab(expression('log'[2]*'(Expression)')) +
    ylab("Time-Replicate") +
    scale_fill_manual(values = brewer.pal(5, "YlOrRd"))
}
#Abundance min and max
ab <- dep |>
  map(SummarizedExperiment::assay) |>
  map(as.matrix) |>
  map(as.vector)

ab |>
  map(min) #10.45

ab |>
  map(max) #24.4

## Violin plots
data_get_results <- results_anno |> 
  map(select,name,ID,name_anno,
      T2_vs_T0_p.val,T4_vs_T0_p.val,      
      T6_vs_T0_p.val,T8_vs_T0_p.val,      
      T2_vs_T0_p.adj,T4_vs_T0_p.adj,      
      T6_vs_T0_p.adj,T8_vs_T0_p.adj,      
      T2_vs_T0_significant,T4_vs_T0_significant,
      T6_vs_T0_significant,T8_vs_T0_significant,
      significant,T2_vs_T0_ratio,   
      T4_vs_T0_ratio,T6_vs_T0_ratio,      
      T8_vs_T0_ratio,T0_centered,         
      T2_centered,T4_centered,         
      T6_centered,T8_centered) |>   
  map(select,-starts_with("T4_vs_T2")) |> 
  map(select,-starts_with("T6_vs_T4")) |> 
  map(select,-starts_with("T8_vs_T6")) |> 
  map(pivot_longer,
      names_to = "category", values_to = "values", c(4:15, 17:20)) |> 
  map(separate, category, sep = "_", into = c("C1", "C2", "C3", "type")) |> 
  map(mutate, Time = paste0(C1, "_", C2, "_", C3)) |>
  map(select, -c("C1", "C2", "C3", "ID")) |> 
  map(rename, Sig.Any = significant) |> 
  map(pivot_wider, names_from = type, values_from = values) |>
  map(mutate, neg.log.p.val = -log(p.val)) |>
  map(mutate, Time = ifelse(Time == "T2_vs_T0", "T2 vs T0",
                            ifelse(Time == "T4_vs_T0", "T4 vs T0",
                                   ifelse(Time == "T6_vs_T0", "T6 vs T0",
                                          ifelse(Time == "T8_vs_T0", "T8 vs T0","NA")))))

violin_plots <- list()
for (i in seq_along(data_get_results)) {
  violin_plots[[i]] <- data_get_results[[i]] |>
    ggplot(aes(x = ratio, y = Time,
               fill = Time)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = "white") +
    theme_minimal() +
    scale_y_discrete(limits = rev) +
    theme(legend.position = "none") +
    ylab("Contrast") +
    xlab(expression('log'[2]*'(Expression T'[X]*'/Expression T'[0]*')')) +
    scale_fill_manual(breaks = c("T2 vs T0", "T4 vs T0", "T6 vs T0", "T8 vs T0"),
                      values = c(brewer.pal(5, "YlOrRd")[2:5]))
}

#Figure 3.
for (p in seq_along(plot_pcas)){
  fig3 <- wrap_plots(distribution_plots[[1]],violin_plots[[1]])+
    plot_annotation(tag_levels = "A", tag_suffix = ".") &
    theme(plot.tag.position = c(0, 1),
          plot.tag = element_text(size = 15, hjust = 0, vjust = 0))
  
  ggsave(
    filename = paste0("proteomics/img_IndependentSets/indset.fig3.",file_set[p],".pdf"),
    plot = fig3,
    width = 9,
    height = 7
  )
  
  ggsave(
    filename = paste0("proteomics/img_IndependentSets/indset.fig3.",file_set[p],".png"),
    plot = fig3,
    width = 9,
    height = 7,dpi=500
  )
  
}  

data_stats <- data_get_results |>
  map(group_by, Time) |>
  map(summarise, min.rat = min(ratio),
      max.rat = max(ratio),
      sd.rat = sd(ratio),
      mean.rat = mean(ratio),
      cv.rat = sd(ratio) / mean(ratio))

df_r <- data_get_results|>
  map(mutate, fold_change = 2^ratio)


# What are the most abundant proteins?
for (p in seq_along(results_anno)){
  
mostabundant <- results_anno[[p]] |>
  select(name,name_anno,T0_1,T0_2,T0_3,
         T2_1,T2_2,T2_3,           
         T4_1,T4_2,T4_3,                
         T6_1,T6_2,T6_3,                
         T8_1,T8_2,T8_3) |> 
  pivot_longer(names_to = "TR", values_to = "exp", 3:17) |> 
  group_by(name_anno) |>
  summarise(max.exp = max(exp)) |>
  arrange(desc(max.exp)) |>
  mutate(max.exp=round(max.exp,2))|> 
  head(20)

mostabundant |>
  write_csv(paste0("proteomics/DEP_results_IndependentSets/mostabundant_",file_set[p],".csv"))
}

# P value histograms, per time point.
pvalue_hist <- list()
for (i in seq_along(data_get_results)) {
  pvalue_hist[[i]] <- data_get_results[[i]] |>
    ggplot(aes(x = p.adj, fill = Time, col = Time)) +
    geom_density() +
    #geom_histogram(binwidth = 0.001, position = "identity") +
    #xlim(0, 0.5) +
    #ylim(0, 20) +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_fill_manual(values = c(brewer.pal(5, "YlOrRd")[2:5])) +
    facet_grid(Time~.)+
    xlab("Adjusted p-values")
}

foldchange_hist <- list()
for (i in seq_along(data_get_results)) {
  foldchange_hist[[i]] <- data_get_results[[i]] |>
    ggplot(aes(x = ratio, fill = Time, col = Time)) +
    geom_density() +
    #geom_histogram(binwidth = 0.001, position = "identity") +
    #xlim(0, 0.5) +
    #ylim(0, 20) +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_fill_manual(values = c(brewer.pal(5, "YlOrRd")[2:5])) +
    facet_grid(Time~.)+
    xlab("Fold change")
}

dep |>
  map(plot_cond)

#How many significant results?
dep |>
  map(get_results) |>
  map(filter, significant == TRUE) |>
  map(summarise, n())


# 530, 394, 623, 623
dep_sigonly <- dep |>
  map(get_results) |>
  map(filter, significant == TRUE)

rownames(dep[[3]])
rownames(dep[[4]])

#HEATMAPS.
ht_significant <- dep |>
  map(plot_heatmap, type = "centered", kmeans = TRUE,
      clustering_distance = "pearson",
      k = 4, col_limit = 1, show_row_names = FALSE,
      indicate = c("condition"),
      left_annotation = rowAnnotation(cluster = anno_block(gp = gpar(fill = 1:4),
                                                           labels = c("A", "B", "C","D"))))
#modify color scale
heatmap_list <- ht_significant
names(heatmap_list) <- c("data_mmetsp","data_jgi","data_det","data_sig")

a <- 2
pdf_path <- paste0("proteomics/img_IndependentSets/indset.fig5.",contr,".pdf")
png_path <- paste0("proteomics/img_IndependentSets/indset.fig5.",contr,".png")

pdf(pdf_path,
    width = 3.5*a,
    height = 4*a)
heatmap_list
dev.off()

png(png_path,
    width = 3.5 * a,
    height = 4 * a,
    units = "in",
    res = 300)
heatmap_list[[3]]
dev.off()

# Obtain clusters
heatmap_list |>
  map(draw)
list.genes <- heatmap_list |>
  map(row_order)
# How many genes are significant in each set?
list.genes |>
  map(unlist) |>
  map(sort) |>
  map(unique) |>
  map(length) 
#mmetsp: 530, jgi: 394, det: 623, sig: 623

cluster1 <- list.genes |>
  map(`[`, c("1")) |>
  map(unlist, use.names = FALSE) |>
  map(as.character)
cluster2 <- list.genes |>
  map(`[`, c("2")) |>
  map(unlist, use.names = FALSE) |>
  map(as.character)
cluster3 <- list.genes |>
  map(`[`, c("3")) |>
  map(unlist, use.names = FALSE) |>
  map(as.character)
cluster4 <- list.genes |>
  map(`[`, c("4")) |>
  map(unlist, use.names = FALSE) |>
  map(as.character)

# How many genes in each cluster in the mmetsp?
cluster1 |> map(length) #vs_t0:77, 158, 80, 91
cluster2 |> map(length) #vs_t0:148, 65, 162, 177
cluster3 |> map(length) #vs_t0:135, 82, 146, 160
cluster4 |> map(length) #vs_t0:170, 89, 183, 195

# Volcano plots (black and grey)
a <- dep |>
  map(plot_volcano, contrast = "T2_vs_T0", label_size = 2, add_names = TRUE)
b <- dep |>
  map(plot_volcano, contrast = "T4_vs_T0", label_size = 2, add_names = TRUE)
c <- dep |>
  map(plot_volcano, contrast = "T6_vs_T0", label_size = 2, add_names = TRUE)
d <- dep |>
  map(plot_volcano, contrast = "T8_vs_T0", label_size = 2, add_names = TRUE)
volcanobw_plots <- pmap(list(a, b, c, d), wrap_plots, ncol = 2, nrow = 2)
wrap_plots(a[[1]],b[[1]],
           c[[1]],d[[1]],nrow=2)
wrap_plots(a[[2]],b[[2]],
           c[[2]],d[[2]],nrow=2)
wrap_plots(a[[3]],b[[3]],
           c[[3]],d[[3]],nrow=2)
wrap_plots(a[[4]],b[[4]],
           c[[4]],d[[4]],nrow=2)

# Volcano plots (color by cluster)
data_get_results <- dep |>
  map(get_results)
row_names <- data_get_results|>
  map(rownames) |>
  map(as.data.frame) |>
  map(rename, rownumber = ".x[[i]]")
data_get_results <- map2(data_get_results, row_names, cbind)

genes <- list()
data_get_results1 <- list()
for (i in seq_along(data_get_results)) {
  genes[[i]] <- data_get_results[[i]] |>
    mutate(cluster = ifelse(rownumber %in% cluster1[[i]], "Cluster 1",
                            ifelse(rownumber %in% cluster2[[i]], "Cluster 2",
                                   ifelse(rownumber %in% cluster3[[i]], "Cluster 3",
                                          ifelse(rownumber %in% cluster4[[i]], "Cluster 4","not significant"))))) |>
    select(name, cluster)
  
  data_get_results1[[i]] <- data_get_results[[i]] |>
    mutate(cluster = ifelse(rownumber %in% cluster1[[i]], "Cluster 1",
                            ifelse(rownumber %in% cluster2[[i]], "Cluster 2",
                                   ifelse(rownumber %in% cluster3[[i]], "Cluster 3",
                                          ifelse(rownumber %in% cluster4[[i]], "Cluster 4","not significant")))))
}

data_get_results2 <- data_get_results1 |> 
  map(select,-starts_with("T4_vs_T2")) |> 
  map(select,-starts_with("T6_vs_T4")) |> 
  map(select,-starts_with("T8_vs_T6")) |> 
  map(pivot_longer,names_to = "category", values_to = "values", c(3:14, 16:19)) |>
  map(separate, category, sep = "_", into = c("C1", "C2", "C3", "type")) |>
  map(mutate, Time = paste0(C1, "_", C2, "_", C3)) |>
  map(select, -c("C1", "C2", "C3", "rownumber", "ID")) |>
  map(rename,any.significant=significant) |> 
  map(pivot_wider, names_from = type, values_from = values) |>
  map(mutate, neg.log.p.val = -log(p.val,10)) |>
  map(mutate, Time = ifelse(Time == "T2_vs_T0", "T2 vs T0",
                            ifelse(Time == "T4_vs_T0", "T4 vs T0",
                                   ifelse(Time == "T6_vs_T0", "T6 vs T0",
                                          ifelse(Time == "T8_vs_T0", "T8 vs T0", "NA")))))

data_sig_false <- data_get_results2|>
  map(filter, any.significant == FALSE)
data_sig_true <- data_get_results2|>
  map(filter, any.significant == TRUE)
thr <- data_sig_false |>
  map(group_by, Time) |>
  map(summarise, max.nlp = max(neg.log.p.val))

volcano_plots <- list()
for (i in seq_along(data_sig_false)) {
  volcano_plots[[i]] <- ggplot() +
    geom_point(aes(x = ratio, y = neg.log.p.val,
                   col = cluster),
               data = data_sig_false[[i]]) +
    geom_point(aes(x = ratio, y = neg.log.p.val,
                   col = cluster),
               data = data_sig_true[[i]]) +
    facet_wrap(.~Time) +
    geom_hline(data = thr[[i]], aes(yintercept = max.nlp),linetype='dotted') +
    geom_vline(xintercept=0)+
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_color_manual(breaks = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "not significant"),
                       values = c(brewer.pal(5, "Set1"), "grey"))+
    xlab(expression('log'[2]*'(Fold change)'))+
    ylab(expression('-log'[10]*'(p value)'))+
    xlim(-1.5, +1.5)+
    ylim(-1,16)+
    geom_vline(xintercept=log2(fold_threshold),linetype='dotted')+
    geom_vline(xintercept=log2(inverse_fold),linetype='dotted')
  
}

a <- 2
pdf(paste0("proteomics/img_IndependentSets/indset.fig6_",contr,".pdf"),
    width = 4*a,
    height = 4*a)
volcano_plots #Figure 6
dev.off()

png_path <- paste0("proteomics/img_IndependentSets/indset.fig6.",contr,".png")
png(png_path,
    width = 3.5 * a,
    height = 4 * a,
    units = "in",
    res = 300)
volcano_plots[[3]]
dev.off()

# Merge cluster to annotations.
cluster_anno <- data_get_results1 |> 
  map(select,name,cluster)

results_anno <- map2(results_anno,cluster_anno,left_join,by="name",suffix = c("", "_n")) |> 
  map(select,!ends_with("_n"))

results_anno_onlysig <- results_anno |>
  map(filter, significant == TRUE)
results_anno_onlysig <- map2(results_anno_onlysig,cluster_anno,
                             left_join,by="name",suffix = c("", "_n")) |> 
  map(select,!ends_with("_n"))

results_anno_3s <- results_anno
results_anno_onlysig_3s <- results_anno_onlysig
# 
# # Protein dynamics of MMETSP set
# len <- dim(dep[[1]])[1]
# gr <- round(len / 200)
# list_dynamics <- list()
# seq1 <- seq(1, len, by = gr)
# seq2 <- seq1 + 16
# 
# # Protein dynamics of all proteins.
# for (j in seq_along(seq1)) {
#   list_dynamics[[j]] <- plot_single(dep[[1]],
#                                     proteins = rownames(dep[[1]])[seq1[j]:seq2[j]],
#                                     type = "centered") +
#     theme_minimal()
# }
# 
# pdf("proteomics/img_IndependentSets/proteindynamics_mmetsp.pdf")
# list_dynamics
# dev.off()
# 
# # Protein dynamics of JGI set
# len <- dim(dep[[2]])[1]
# gr <- round(len / 200)
# list_dynamics <- list()
# seq1 <- seq(1, len, by = gr)
# seq2 <- seq1 + 16
# 
# # Protein dynamics of all proteins.
# for (j in seq_along(seq1)) {
#   list_dynamics[[j]] <- plot_single(dep[[2]],
#                                     proteins = rownames(dep[[2]])[seq1[j]:seq2[j]],
#                                     type = "centered") +
#     theme_minimal()
# }
# 
# pdf("proteomics/img_IndependentSets/proteindynamics_jgi.pdf")
# list_dynamics
# dev.off()
# 
# 
# # Protein dynamics of ALLDET set
# len <- dim(dep[[3]])[1]
# gr <- round(len / 200)
# list_dynamics <- list()
# seq1 <- seq(1, len, by = gr)
# seq2 <- seq1 + 16
# 
# # Protein dynamics of all proteins.
# for (j in seq_along(seq1)) {
#   list_dynamics[[j]] <- plot_single(dep[[3]],
#                                     proteins = rownames(dep[[3]])[seq1[j]:seq2[j]],
#                                     type = "centered") +
#     theme_minimal()
# }
# 
# pdf("proteomics/img_IndependentSets/proteindynamics_allsig.pdf")
# list_dynamics
# dev.off()


## Save tables, with annotations.
df_wide <- dep |>
  map(get_df_wide)
df_long <- dep |>
  map(get_df_long)

df_wide <- map2(df_wide, genes, left_join, by = "name")
df_long <- map2(df_long, genes, left_join, by = "name")

for (i in seq_along(df_wide)) {
  df_wide[[i]] |>
    write_csv(file.path("proteomics/DEP_results_IndependentSets", paste0("df_wide.",file_set[i],".",contr,".csv")))
  df_long[[i]] |>
    write_csv(file.path("proteomics/DEP_results_IndependentSets", paste0("df_long.", file_set[i],".",contr,".csv")))
}

# With annotations.
results_anno
results_anno_onlysig
results_anno_3s
results_anno_onlysig_3s
for (i in seq_along(df_wide)) {
  results_anno[[i]] |>
    write_csv(file.path("proteomics/DEP_results_IndependentSets", paste0("df_anno.", file_set[i],".",contr,".csv")))
  results_anno_onlysig[[i]] |>
    write_csv(file.path("proteomics/DEP_results_IndependentSets", paste0("df_anno_filtsig.", file_set[i],".",contr,".csv")))
}

#Save files with annotations.
if (contr == "t0") {
  threesets_concat <- rbind(results_anno[[1]], results_anno[[2]]) |>
    group_by(name_anno) |>
    summarise(
      "RAW.T0-1" =  paste0(`RAW.T0-1`, collapse = "; "),
      "RAW.T0-2" =  paste0(`RAW.T0-2`, collapse = "; "),
      "RAW.T0-3" =  paste0(`RAW.T0-3`, collapse = "; "),
      "RAW.T2-1" =  paste0(`RAW.T2-1`, collapse = "; "),
      "RAW.T2-2" =  paste0(`RAW.T2-2`, collapse = "; "),
      "RAW.T2-3" =  paste0(`RAW.T2-3`, collapse = "; "),
      "RAW.T4-1" =  paste0(`RAW.T4-1`, collapse = "; "),
      "RAW.T4-2" =  paste0(`RAW.T4-2`, collapse = "; "),
      "RAW.T4-3" =  paste0(`RAW.T4-3`, collapse = "; "),
      "RAW.T6-1" =  paste0(`RAW.T6-1`, collapse = "; "),
      "RAW.T6-2" =  paste0(`RAW.T6-2`, collapse = "; "),
      "RAW.T6-3" =  paste0(`RAW.T6-3`, collapse = "; "),
      "RAW.T8-1" =  paste0(`RAW.T8-1`, collapse = "; "),
      "RAW.T8-2" =  paste0(`RAW.T8-2`, collapse = "; "),
      "RAW.T8-3" =  paste0(`RAW.T8-3`, collapse = "; "),
      
      "BATCH_CORRECT.T0-1" =  paste0(`BATCH_CORRECT.T0-1`, collapse = "; "),
      "BATCH_CORRECT.T0-2" =  paste0(`BATCH_CORRECT.T0-2`, collapse = "; "),
      "BATCH_CORRECT.T0-3" =  paste0(`BATCH_CORRECT.T0-3`, collapse = "; "),
      "BATCH_CORRECT.T2-1" =  paste0(`BATCH_CORRECT.T2-1`, collapse = "; "),
      "BATCH_CORRECT.T2-2" =  paste0(`BATCH_CORRECT.T2-2`, collapse = "; "),
      "BATCH_CORRECT.T2-3" =  paste0(`BATCH_CORRECT.T2-3`, collapse = "; "),
      "BATCH_CORRECT.T4-1" =  paste0(`BATCH_CORRECT.T4-1`, collapse = "; "),
      "BATCH_CORRECT.T4-2" =  paste0(`BATCH_CORRECT.T4-2`, collapse = "; "),
      "BATCH_CORRECT.T4-3" =  paste0(`BATCH_CORRECT.T4-3`, collapse = "; "),
      "BATCH_CORRECT.T6-1" =  paste0(`BATCH_CORRECT.T6-1`, collapse = "; "),
      "BATCH_CORRECT.T6-2" =  paste0(`BATCH_CORRECT.T6-2`, collapse = "; "),
      "BATCH_CORRECT.T6-3" =  paste0(`BATCH_CORRECT.T6-3`, collapse = "; "),
      "BATCH_CORRECT.T8-1" =  paste0(`BATCH_CORRECT.T8-1`, collapse = "; "),
      "BATCH_CORRECT.T8-2" =  paste0(`BATCH_CORRECT.T8-2`, collapse = "; "),
      "BATCH_CORRECT.T8-3" =  paste0(`BATCH_CORRECT.T8-3`, collapse = "; "),
      
      "T0_1" = paste0(T0_1, collapse = "; "),
      "T0_2" = paste0(T0_2, collapse = "; "),
      "T0_3" = paste0(T0_3, collapse = "; "),
      "T2_1" = paste0(T2_1, collapse = "; "),
      "T2_2" = paste0(T2_2, collapse = "; "),
      "T2_3" = paste0(T2_3, collapse = "; "),
      "T4_1" = paste0(T4_1, collapse = "; "),
      "T4_2" = paste0(T4_2, collapse = "; "),
      "T4_3" = paste0(T4_3, collapse = "; "),
      "T6_1" = paste0(T6_1, collapse = "; "),
      "T6_2" = paste0(T6_2, collapse = "; "),
      "T6_3" = paste0(T6_3, collapse = "; "),
      "T8_1" = paste0(T8_1, collapse = "; "),
      "T8_2" = paste0(T8_2, collapse = "; "),
      "T8_3" = paste0(T8_3, collapse = "; "),
      
      "T2_vs_T0_ratio" = paste0(T2_vs_T0_ratio, collapse = "; "),
      "T4_vs_T0_ratio" = paste0(T4_vs_T0_ratio, collapse = "; "),
      "T6_vs_T0_ratio" = paste0(T6_vs_T0_ratio, collapse = "; "),
      "T8_vs_T0_ratio" = paste0(T8_vs_T0_ratio, collapse = "; "),
      
      "T2_centered" = paste0(T2_centered, collapse = "; "),
      "T4_centered" = paste0(T4_centered, collapse = "; "),
      "T6_centered" = paste0(T6_centered, collapse = "; "),
      "T8_centered" = paste0(T8_centered, collapse = "; "),
      
      "T2_vs_T0_p.val" = paste0(T2_vs_T0_p.val, collapse = "; "),
      "T4_vs_T0_p.val" = paste0(T4_vs_T0_p.val, collapse = "; "),
      "T6_vs_T0_p.val" = paste0(T6_vs_T0_p.val, collapse = "; "),
      "T8_vs_T0_p.val" = paste0(T8_vs_T0_p.val, collapse = "; "),
      
      "T2_vs_T0_p.adj" = paste0(T2_vs_T0_p.adj, collapse = "; "),
      "T4_vs_T0_p.adj" = paste0(T4_vs_T0_p.adj, collapse = "; "),
      "T6_vs_T0_p.adj" = paste0(T6_vs_T0_p.adj, collapse = "; "),
      "T8_vs_T0_p.adj" = paste0(T8_vs_T0_p.adj, collapse = "; "),
      
      "T2_vs_T0_significant" = paste0(T2_vs_T0_significant, collapse = "; "),
      "T4_vs_T0_significant" = paste0(T4_vs_T0_significant, collapse = "; "),
      "T6_vs_T0_significant" = paste0(T6_vs_T0_significant, collapse = "; "),
      "T8_vs_T0_significant" = paste0(T8_vs_T0_significant, collapse = "; "),
      "significant" = paste0(significant, collapse = "; "),
      
      "T2_vs_T0_CI.L" = paste0(T2_vs_T0_CI.L, collapse = "; "),
      "T2_vs_T0_CI.R" = paste0(T2_vs_T0_CI.R, collapse = "; "),
      "T2_vs_T0_diff" = paste0(T2_vs_T0_diff, collapse = "; "),
      
      "T4_vs_T0_CI.L" = paste0(T4_vs_T0_CI.L, collapse = "; "),
      "T4_vs_T0_CI.R" = paste0(T4_vs_T0_CI.R, collapse = "; "),
      "T4_vs_T0_diff" = paste0(T4_vs_T0_diff, collapse = "; "),
      
      "T6_vs_T0_CI.L" = paste0(T6_vs_T0_CI.L, collapse = "; "),
      "T6_vs_T0_CI.R" = paste0(T6_vs_T0_CI.R, collapse = "; "),
      "T6_vs_T0_diff" = paste0(T6_vs_T0_diff, collapse = "; "),
      
      "T8_vs_T0_CI.L" = paste0(T8_vs_T0_CI.L, collapse = "; "),
      "T8_vs_T0_CI.R" = paste0(T8_vs_T0_CI.R, collapse = "; "),
      "T8_vs_T0_diff" = paste0(T8_vs_T0_diff, collapse = "; "),
      
      "FC_T2_vs_T0" = paste0(FC_T2_vs_T0, collapse = "; "),
      "FC_T4_vs_T0" = paste0(FC_T4_vs_T0, collapse = "; "),
      "FC_T6_vs_T0" = paste0(FC_T6_vs_T0, collapse = "; "),
      "FC_T8_vs_T0" = paste0(FC_T8_vs_T0, collapse = "; "),
      
      "T2v0.reg" = paste0(T2v0.reg, collapse = "; "),
      "T4v0.reg" = paste0(T4v0.reg, collapse = "; "),
      "T6v0.reg" = paste0(T6v0.reg, collapse = "; "),
      "T8v0.reg" = paste0(T8v0.reg, collapse = "; "),
      
      "orthocluster" = paste0(orthocluster, collapse = "; "),
      "Aplke1" = paste0(Aplke1, collapse = "; "),
      "Aurli1" = paste0(Aurli1, collapse = "; "),
      "Honfer1" = paste0(Honfer1, collapse = "; "),
      "LabyHa" = paste0(LabyHa, collapse = "; "),
      "Schag1" = paste0(Schag1, collapse = "; "),
      
      "num.species" = paste0(num.species, collapse = "; "),
      "conserved" = paste0(conserved, collapse = "; "),
      
      "above_median_exp_t0" = paste0(above_median_exp_t0, collapse = "; "),
      "Annotations" = paste0(Annotations, collapse = "; "),
      "genbank" = paste0(genbank, collapse = "; "),
      "Uniprot_Hon" = paste0(Uniprot_Hon, collapse = "; "),
      "gotermId" = paste0(gotermId, collapse = "; "),
      "goName" = paste0(goName, collapse = "; "),
      "gotermType" = paste0(gotermType, collapse = "; "),
      "goAcc" = paste0(goAcc, collapse = "; "),
      "kogdefline" = paste0(kogdefline, collapse = "; "),
      "kogClass" = paste0(kogClass, collapse = "; "),
      "GOcompartment" = paste0(GOcompartment, collapse = "; "),
      "Sig_Pfam_info" = paste0(Sig_Pfam_info, collapse = "; "),
      "Sig_Pfam_number" = paste0(Sig_Pfam_number, collapse = "; "),
      "Pfam_info" = paste0(Pfam_info, collapse = "; "),
      "Pfam_IPS" = paste0(Pfam_IPS, collapse = "; "),
      "PANTHER_info" = paste0(PANTHER_info, collapse = "; "),
      "PANTHER_IPS" = paste0(PANTHER_IPS, collapse = "; "),
      "multiloc1" = paste0(multiloc1, collapse = "; "),
      "multiloc2" = paste0(multiloc2, collapse = "; "),
      "multiloc3" = paste0(multiloc3, collapse = "; "),
      "KOG_number" = paste0(KOG_number, collapse = "; "),
      "KOG_info" = paste0(KOG_info, collapse = "; "),
      "Keg" = paste0(Keg, collapse = "; "),
      "ko_values" = paste0(ko_values, collapse = "; "),
      "seq_name" = paste0(seq_name, collapse = "; "),
      "UP_AND468" = paste0(UP_AND468, collapse = "; "),
      "UP_OR468" = paste0(UP_OR468, collapse = "; "),
      "DOWN_AND468" = paste0(DOWN_AND468, collapse = "; "),
      "DOWN_OR468" = paste0(DOWN_OR468, collapse = "; "),
      "Increasing" = paste0(Increasing, collapse = "; "),
      "table" = paste0(table, collapse = "; "),
      "set" = paste0(set, collapse = "; "))
} else if (contr == "tx") {
  
  threesets_concat <- rbind(results_anno[[1]], results_anno[[2]]) |>
    group_by(name_anno) |>
    summarise(
      "RAW.T0-1" =  paste0(`RAW.T0-1`, collapse = "; "),
      "RAW.T0-2" =  paste0(`RAW.T0-2`, collapse = "; "),
      "RAW.T0-3" =  paste0(`RAW.T0-3`, collapse = "; "),
      "RAW.T2-1" =  paste0(`RAW.T2-1`, collapse = "; "),
      "RAW.T2-2" =  paste0(`RAW.T2-2`, collapse = "; "),
      "RAW.T2-3" =  paste0(`RAW.T2-3`, collapse = "; "),
      "RAW.T4-1" =  paste0(`RAW.T4-1`, collapse = "; "),
      "RAW.T4-2" =  paste0(`RAW.T4-2`, collapse = "; "),
      "RAW.T4-3" =  paste0(`RAW.T4-3`, collapse = "; "),
      "RAW.T6-1" =  paste0(`RAW.T6-1`, collapse = "; "),
      "RAW.T6-2" =  paste0(`RAW.T6-2`, collapse = "; "),
      "RAW.T6-3" =  paste0(`RAW.T6-3`, collapse = "; "),
      "RAW.T8-1" =  paste0(`RAW.T8-1`, collapse = "; "),
      "RAW.T8-2" =  paste0(`RAW.T8-2`, collapse = "; "),
      "RAW.T8-3" =  paste0(`RAW.T8-3`, collapse = "; "),
      
      "BATCH_CORRECT.T0-1" =  paste0(`BATCH_CORRECT.T0-1`, collapse = "; "),
      "BATCH_CORRECT.T0-2" =  paste0(`BATCH_CORRECT.T0-2`, collapse = "; "),
      "BATCH_CORRECT.T0-3" =  paste0(`BATCH_CORRECT.T0-3`, collapse = "; "),
      "BATCH_CORRECT.T2-1" =  paste0(`BATCH_CORRECT.T2-1`, collapse = "; "),
      "BATCH_CORRECT.T2-2" =  paste0(`BATCH_CORRECT.T2-2`, collapse = "; "),
      "BATCH_CORRECT.T2-3" =  paste0(`BATCH_CORRECT.T2-3`, collapse = "; "),
      "BATCH_CORRECT.T4-1" =  paste0(`BATCH_CORRECT.T4-1`, collapse = "; "),
      "BATCH_CORRECT.T4-2" =  paste0(`BATCH_CORRECT.T4-2`, collapse = "; "),
      "BATCH_CORRECT.T4-3" =  paste0(`BATCH_CORRECT.T4-3`, collapse = "; "),
      "BATCH_CORRECT.T6-1" =  paste0(`BATCH_CORRECT.T6-1`, collapse = "; "),
      "BATCH_CORRECT.T6-2" =  paste0(`BATCH_CORRECT.T6-2`, collapse = "; "),
      "BATCH_CORRECT.T6-3" =  paste0(`BATCH_CORRECT.T6-3`, collapse = "; "),
      "BATCH_CORRECT.T8-1" =  paste0(`BATCH_CORRECT.T8-1`, collapse = "; "),
      "BATCH_CORRECT.T8-2" =  paste0(`BATCH_CORRECT.T8-2`, collapse = "; "),
      "BATCH_CORRECT.T8-3" =  paste0(`BATCH_CORRECT.T8-3`, collapse = "; "),
      
      "T0_1" = paste0(T0_1, collapse = "; "),
      "T0_2" = paste0(T0_2, collapse = "; "),
      "T0_3" = paste0(T0_3, collapse = "; "),
      "T2_1" = paste0(T2_1, collapse = "; "),
      "T2_2" = paste0(T2_2, collapse = "; "),
      "T2_3" = paste0(T2_3, collapse = "; "),
      "T4_1" = paste0(T4_1, collapse = "; "),
      "T4_2" = paste0(T4_2, collapse = "; "),
      "T4_3" = paste0(T4_3, collapse = "; "),
      "T6_1" = paste0(T6_1, collapse = "; "),
      "T6_2" = paste0(T6_2, collapse = "; "),
      "T6_3" = paste0(T6_3, collapse = "; "),
      "T8_1" = paste0(T8_1, collapse = "; "),
      "T8_2" = paste0(T8_2, collapse = "; "),
      "T8_3" = paste0(T8_3, collapse = "; "),
      
      "T2_vs_T0_ratio" = paste0(T2_vs_T0_ratio, collapse = "; "),
      "T4_vs_T0_ratio" = paste0(T4_vs_T0_ratio, collapse = "; "),
      "T6_vs_T0_ratio" = paste0(T6_vs_T0_ratio, collapse = "; "),
      "T8_vs_T0_ratio" = paste0(T8_vs_T0_ratio, collapse = "; "),
      "T4_vs_T2_ratio" = paste0(T4_vs_T2_ratio, collapse = "; "),
      "T6_vs_T4_ratio" = paste0(T6_vs_T4_ratio, collapse = "; "),
      "T8_vs_T6_ratio" = paste0(T8_vs_T6_ratio, collapse = "; "),
      
      "T2_centered" = paste0(T2_centered, collapse = "; "),
      "T4_centered" = paste0(T4_centered, collapse = "; "),
      "T6_centered" = paste0(T6_centered, collapse = "; "),
      "T8_centered" = paste0(T8_centered, collapse = "; "),
      
      "T2_vs_T0_p.val" = paste0(T2_vs_T0_p.val, collapse = "; "),
      "T4_vs_T0_p.val" = paste0(T4_vs_T0_p.val, collapse = "; "),
      "T6_vs_T0_p.val" = paste0(T6_vs_T0_p.val, collapse = "; "),
      "T8_vs_T0_p.val" = paste0(T8_vs_T0_p.val, collapse = "; "),
      "T4_vs_T2_p.val" = paste0(T4_vs_T2_p.val, collapse = "; "),
      "T6_vs_T4_p.val" = paste0(T6_vs_T4_p.val, collapse = "; "),
      "T8_vs_T6_p.val" = paste0(T8_vs_T6_p.val, collapse = "; "),
      
      "T2_vs_T0_p.adj" = paste0(T2_vs_T0_p.adj, collapse = "; "),
      "T4_vs_T0_p.adj" = paste0(T4_vs_T0_p.adj, collapse = "; "),
      "T6_vs_T0_p.adj" = paste0(T6_vs_T0_p.adj, collapse = "; "),
      "T8_vs_T0_p.adj" = paste0(T8_vs_T0_p.adj, collapse = "; "),
      "T4_vs_T2_p.adj" = paste0(T4_vs_T2_p.adj, collapse = "; "),
      "T6_vs_T4_p.adj" = paste0(T6_vs_T4_p.adj, collapse = "; "),
      "T8_vs_T6_p.adj" = paste0(T8_vs_T6_p.adj, collapse = "; "),
      
      "T2_vs_T0_significant" = paste0(T2_vs_T0_significant, collapse = "; "),
      "T4_vs_T0_significant" = paste0(T4_vs_T0_significant, collapse = "; "),
      "T6_vs_T0_significant" = paste0(T6_vs_T0_significant, collapse = "; "),
      "T8_vs_T0_significant" = paste0(T8_vs_T0_significant, collapse = "; "),
      "T4_vs_T2_significant" = paste0(T4_vs_T2_significant, collapse = "; "),
      "T6_vs_T4_significant" = paste0(T6_vs_T4_significant, collapse = "; "),
      "T8_vs_T6_significant" = paste0(T8_vs_T6_significant, collapse = "; "),
      "significant" = paste0(significant, collapse = "; "),
      
      "T2_vs_T0_CI.L" = paste0(T2_vs_T0_CI.L, collapse = "; "),
      "T2_vs_T0_CI.R" = paste0(T2_vs_T0_CI.R, collapse = "; "),
      "T2_vs_T0_diff" = paste0(T2_vs_T0_diff, collapse = "; "),
      
      "T4_vs_T0_CI.L" = paste0(T4_vs_T0_CI.L, collapse = "; "),
      "T4_vs_T0_CI.R" = paste0(T4_vs_T0_CI.R, collapse = "; "),
      "T4_vs_T0_diff" = paste0(T4_vs_T0_diff, collapse = "; "),
      
      "T6_vs_T0_CI.L" = paste0(T6_vs_T0_CI.L, collapse = "; "),
      "T6_vs_T0_CI.R" = paste0(T6_vs_T0_CI.R, collapse = "; "),
      "T6_vs_T0_diff" = paste0(T6_vs_T0_diff, collapse = "; "),
      
      "T8_vs_T0_CI.L" = paste0(T8_vs_T0_CI.L, collapse = "; "),
      "T8_vs_T0_CI.R" = paste0(T8_vs_T0_CI.R, collapse = "; "),
      "T8_vs_T0_diff" = paste0(T8_vs_T0_diff, collapse = "; "),
      
      "T4_vs_T2_CI.L" = paste0(T4_vs_T2_CI.L, collapse = "; "),
      "T4_vs_T2_CI.R" = paste0(T4_vs_T2_CI.R, collapse = "; "),
      "T4_vs_T2_diff" = paste0(T4_vs_T2_diff, collapse = "; "),
      
      "T6_vs_T4_CI.L" = paste0(T6_vs_T4_CI.L, collapse = "; "),
      "T6_vs_T4_CI.R" = paste0(T6_vs_T4_CI.R, collapse = "; "),
      "T6_vs_T4_diff" = paste0(T6_vs_T4_diff, collapse = "; "),
      
      "T8_vs_T6_CI.L" = paste0(T8_vs_T6_CI.L, collapse = "; "),
      "T8_vs_T6_CI.R" = paste0(T8_vs_T6_CI.R, collapse = "; "),
      "T8_vs_T6_diff" = paste0(T8_vs_T6_diff, collapse = "; "),
      
      "FC_T2_vs_T0" = paste0(FC_T2_vs_T0, collapse = "; "),
      "FC_T4_vs_T0" = paste0(FC_T4_vs_T0, collapse = "; "),
      "FC_T6_vs_T0" = paste0(FC_T6_vs_T0, collapse = "; "),
      "FC_T8_vs_T0" = paste0(FC_T8_vs_T0, collapse = "; "),
      "FC_T4_vs_T2" = paste0(FC_T4_vs_T2, collapse = "; "),
      "FC_T6_vs_T4" = paste0(FC_T6_vs_T4, collapse = "; "),
      "FC_T8_vs_T6" = paste0(FC_T8_vs_T6, collapse = "; "),
      
      "T2v0.reg" = paste0(T2v0.reg, collapse = "; "),
      "T4v0.reg" = paste0(T4v0.reg, collapse = "; "),
      "T6v0.reg" = paste0(T6v0.reg, collapse = "; "),
      "T8v0.reg" = paste0(T8v0.reg, collapse = "; "),
      "T4v2.reg" = paste0(T4v2.reg, collapse = "; "),
      "T6v4.reg" = paste0(T6v4.reg, collapse = "; "),
      "T8v6.reg" = paste0(T8v6.reg, collapse = "; "),
      
      "orthocluster" = paste0(orthocluster, collapse = "; "),
      "Aplke1" = paste0(Aplke1, collapse = "; "),
      "Aurli1" = paste0(Aurli1, collapse = "; "),
      "Honfer1" = paste0(Honfer1, collapse = "; "),
      "LabyHa" = paste0(LabyHa, collapse = "; "),
      "Schag1" = paste0(Schag1, collapse = "; "),
      
      "num.species" = paste0(num.species, collapse = "; "),
      "conserved" = paste0(conserved, collapse = "; "),
      
      "above_median_exp_t0" = paste0(above_median_exp_t0, collapse = "; "),
      "Annotations" = paste0(Annotations, collapse = "; "),
      "genbank" = paste0(genbank, collapse = "; "),
      "Uniprot_Hon" = paste0(Uniprot_Hon, collapse = "; "),
      "gotermId" = paste0(gotermId, collapse = "; "),
      "goName" = paste0(goName, collapse = "; "),
      "gotermType" = paste0(gotermType, collapse = "; "),
      "goAcc" = paste0(goAcc, collapse = "; "),
      "kogdefline" = paste0(kogdefline, collapse = "; "),
      "kogClass" = paste0(kogClass, collapse = "; "),
      "GOcompartment" = paste0(GOcompartment, collapse = "; "),
      "Sig_Pfam_info" = paste0(Sig_Pfam_info, collapse = "; "),
      "Sig_Pfam_number" = paste0(Sig_Pfam_number, collapse = "; "),
      "Pfam_info" = paste0(Pfam_info, collapse = "; "),
      "Pfam_IPS" = paste0(Pfam_IPS, collapse = "; "),
      "PANTHER_info" = paste0(PANTHER_info, collapse = "; "),
      "PANTHER_IPS" = paste0(PANTHER_IPS, collapse = "; "),
      "multiloc1" = paste0(multiloc1, collapse = "; "),
      "multiloc2" = paste0(multiloc2, collapse = "; "),
      "multiloc3" = paste0(multiloc3, collapse = "; "),
      "KOG_number" = paste0(KOG_number, collapse = "; "),
      "KOG_info" = paste0(KOG_info, collapse = "; "),
      "Keg" = paste0(Keg, collapse = "; "),
      "ko_values" = paste0(ko_values, collapse = "; "),
      "seq_name" = paste0(seq_name, collapse = "; "),
      "UP_AND468" = paste0(UP_AND468, collapse = "; "),
      "UP_OR468" = paste0(UP_OR468, collapse = "; "),
      "DOWN_AND468" = paste0(DOWN_AND468, collapse = "; "),
      "DOWN_OR468" = paste0(DOWN_OR468, collapse = "; "),
      "UP_T6T4_OR_T4T2" = paste0(UP_T6T4_OR_T4T2, collapse = "; "),
      "Increasing" = paste0(Increasing, collapse = "; "),
      "table" = paste0(table, collapse = "; "),
      "set" = paste0(set, collapse = "; "))
} else {
  # Code to execute if 'set' is neither "t0" nor "tx"
  # Add your code here
  print("Error")
}

## Results (filtering only sig)
results_anno_onlysig <- results_anno |>
  map(filter, significant == TRUE)

if (contr == "t0") {
  threesets_concat_sig <- rbind(results_anno_onlysig[[1]],
                                results_anno_onlysig[[2]]) |>
    group_by(name_anno) |>
    summarise(
      "RAW.T0-1" =  paste0(`RAW.T0-1`, collapse = "; "),
      "RAW.T0-2" =  paste0(`RAW.T0-2`, collapse = "; "),
      "RAW.T0-3" =  paste0(`RAW.T0-3`, collapse = "; "),
      "RAW.T2-1" =  paste0(`RAW.T2-1`, collapse = "; "),
      "RAW.T2-2" =  paste0(`RAW.T2-2`, collapse = "; "),
      "RAW.T2-3" =  paste0(`RAW.T2-3`, collapse = "; "),
      "RAW.T4-1" =  paste0(`RAW.T4-1`, collapse = "; "),
      "RAW.T4-2" =  paste0(`RAW.T4-2`, collapse = "; "),
      "RAW.T4-3" =  paste0(`RAW.T4-3`, collapse = "; "),
      "RAW.T6-1" =  paste0(`RAW.T6-1`, collapse = "; "),
      "RAW.T6-2" =  paste0(`RAW.T6-2`, collapse = "; "),
      "RAW.T6-3" =  paste0(`RAW.T6-3`, collapse = "; "),
      "RAW.T8-1" =  paste0(`RAW.T8-1`, collapse = "; "),
      "RAW.T8-2" =  paste0(`RAW.T8-2`, collapse = "; "),
      "RAW.T8-3" =  paste0(`RAW.T8-3`, collapse = "; "),
      
      "BATCH_CORRECT.T0-1" =  paste0(`BATCH_CORRECT.T0-1`, collapse = "; "),
      "BATCH_CORRECT.T0-2" =  paste0(`BATCH_CORRECT.T0-2`, collapse = "; "),
      "BATCH_CORRECT.T0-3" =  paste0(`BATCH_CORRECT.T0-3`, collapse = "; "),
      "BATCH_CORRECT.T2-1" =  paste0(`BATCH_CORRECT.T2-1`, collapse = "; "),
      "BATCH_CORRECT.T2-2" =  paste0(`BATCH_CORRECT.T2-2`, collapse = "; "),
      "BATCH_CORRECT.T2-3" =  paste0(`BATCH_CORRECT.T2-3`, collapse = "; "),
      "BATCH_CORRECT.T4-1" =  paste0(`BATCH_CORRECT.T4-1`, collapse = "; "),
      "BATCH_CORRECT.T4-2" =  paste0(`BATCH_CORRECT.T4-2`, collapse = "; "),
      "BATCH_CORRECT.T4-3" =  paste0(`BATCH_CORRECT.T4-3`, collapse = "; "),
      "BATCH_CORRECT.T6-1" =  paste0(`BATCH_CORRECT.T6-1`, collapse = "; "),
      "BATCH_CORRECT.T6-2" =  paste0(`BATCH_CORRECT.T6-2`, collapse = "; "),
      "BATCH_CORRECT.T6-3" =  paste0(`BATCH_CORRECT.T6-3`, collapse = "; "),
      "BATCH_CORRECT.T8-1" =  paste0(`BATCH_CORRECT.T8-1`, collapse = "; "),
      "BATCH_CORRECT.T8-2" =  paste0(`BATCH_CORRECT.T8-2`, collapse = "; "),
      "BATCH_CORRECT.T8-3" =  paste0(`BATCH_CORRECT.T8-3`, collapse = "; "),
      
      "T0_1" = paste0(T0_1, collapse = "; "),
      "T0_2" = paste0(T0_2, collapse = "; "),
      "T0_3" = paste0(T0_3, collapse = "; "),
      "T2_1" = paste0(T2_1, collapse = "; "),
      "T2_2" = paste0(T2_2, collapse = "; "),
      "T2_3" = paste0(T2_3, collapse = "; "),
      "T4_1" = paste0(T4_1, collapse = "; "),
      "T4_2" = paste0(T4_2, collapse = "; "),
      "T4_3" = paste0(T4_3, collapse = "; "),
      "T6_1" = paste0(T6_1, collapse = "; "),
      "T6_2" = paste0(T6_2, collapse = "; "),
      "T6_3" = paste0(T6_3, collapse = "; "),
      "T8_1" = paste0(T8_1, collapse = "; "),
      "T8_2" = paste0(T8_2, collapse = "; "),
      "T8_3" = paste0(T8_3, collapse = "; "),
      
      "T2_vs_T0_ratio" = paste0(T2_vs_T0_ratio, collapse = "; "),
      "T4_vs_T0_ratio" = paste0(T4_vs_T0_ratio, collapse = "; "),
      "T6_vs_T0_ratio" = paste0(T6_vs_T0_ratio, collapse = "; "),
      "T8_vs_T0_ratio" = paste0(T8_vs_T0_ratio, collapse = "; "),
      
      "T2_centered" = paste0(T2_centered, collapse = "; "),
      "T4_centered" = paste0(T4_centered, collapse = "; "),
      "T6_centered" = paste0(T6_centered, collapse = "; "),
      "T8_centered" = paste0(T8_centered, collapse = "; "),
      
      "T2_vs_T0_p.val" = paste0(T2_vs_T0_p.val, collapse = "; "),
      "T4_vs_T0_p.val" = paste0(T4_vs_T0_p.val, collapse = "; "),
      "T6_vs_T0_p.val" = paste0(T6_vs_T0_p.val, collapse = "; "),
      "T8_vs_T0_p.val" = paste0(T8_vs_T0_p.val, collapse = "; "),
      
      "T2_vs_T0_p.adj" = paste0(T2_vs_T0_p.adj, collapse = "; "),
      "T4_vs_T0_p.adj" = paste0(T4_vs_T0_p.adj, collapse = "; "),
      "T6_vs_T0_p.adj" = paste0(T6_vs_T0_p.adj, collapse = "; "),
      "T8_vs_T0_p.adj" = paste0(T8_vs_T0_p.adj, collapse = "; "),
      
      "T2_vs_T0_significant" = paste0(T2_vs_T0_significant, collapse = "; "),
      "T4_vs_T0_significant" = paste0(T4_vs_T0_significant, collapse = "; "),
      "T6_vs_T0_significant" = paste0(T6_vs_T0_significant, collapse = "; "),
      "T8_vs_T0_significant" = paste0(T8_vs_T0_significant, collapse = "; "),
      "significant" = paste0(significant, collapse = "; "),
      
      "T2_vs_T0_CI.L" = paste0(T2_vs_T0_CI.L, collapse = "; "),
      "T2_vs_T0_CI.R" = paste0(T2_vs_T0_CI.R, collapse = "; "),
      "T2_vs_T0_diff" = paste0(T2_vs_T0_diff, collapse = "; "),
      
      "T4_vs_T0_CI.L" = paste0(T4_vs_T0_CI.L, collapse = "; "),
      "T4_vs_T0_CI.R" = paste0(T4_vs_T0_CI.R, collapse = "; "),
      "T4_vs_T0_diff" = paste0(T4_vs_T0_diff, collapse = "; "),
      
      "T6_vs_T0_CI.L" = paste0(T6_vs_T0_CI.L, collapse = "; "),
      "T6_vs_T0_CI.R" = paste0(T6_vs_T0_CI.R, collapse = "; "),
      "T6_vs_T0_diff" = paste0(T6_vs_T0_diff, collapse = "; "),
      
      "T8_vs_T0_CI.L" = paste0(T8_vs_T0_CI.L, collapse = "; "),
      "T8_vs_T0_CI.R" = paste0(T8_vs_T0_CI.R, collapse = "; "),
      "T8_vs_T0_diff" = paste0(T8_vs_T0_diff, collapse = "; "),
      
      "FC_T2_vs_T0" = paste0(FC_T2_vs_T0, collapse = "; "),
      "FC_T4_vs_T0" = paste0(FC_T4_vs_T0, collapse = "; "),
      "FC_T6_vs_T0" = paste0(FC_T6_vs_T0, collapse = "; "),
      "FC_T8_vs_T0" = paste0(FC_T8_vs_T0, collapse = "; "),
      
      "T2v0.reg" = paste0(T2v0.reg, collapse = "; "),
      "T4v0.reg" = paste0(T4v0.reg, collapse = "; "),
      "T6v0.reg" = paste0(T6v0.reg, collapse = "; "),
      "T8v0.reg" = paste0(T8v0.reg, collapse = "; "),
      
      "orthocluster" = paste0(orthocluster, collapse = "; "),
      "Aplke1" = paste0(Aplke1, collapse = "; "),
      "Aurli1" = paste0(Aurli1, collapse = "; "),
      "Honfer1" = paste0(Honfer1, collapse = "; "),
      "LabyHa" = paste0(LabyHa, collapse = "; "),
      "Schag1" = paste0(Schag1, collapse = "; "),
      "num.species" = paste0(num.species, collapse = "; "),
      "conserved" = paste0(conserved, collapse = "; "),
      
      "above_median_exp_t0" = paste0(above_median_exp_t0, collapse = "; "),
      "Annotations" = paste0(Annotations, collapse = "; "),
      "genbank" = paste0(genbank, collapse = "; "),
      "Uniprot_Hon" = paste0(Uniprot_Hon, collapse = "; "),
      "gotermId" = paste0(gotermId, collapse = "; "),
      "goName" = paste0(goName, collapse = "; "),
      "gotermType" = paste0(gotermType, collapse = "; "),
      "goAcc" = paste0(goAcc, collapse = "; "),
      "kogdefline" = paste0(kogdefline, collapse = "; "),
      "kogClass" = paste0(kogClass, collapse = "; "),
      "GOcompartment" = paste0(GOcompartment, collapse = "; "),
      "Sig_Pfam_info" = paste0(Sig_Pfam_info, collapse = "; "),
      "Sig_Pfam_number" = paste0(Sig_Pfam_number, collapse = "; "),
      "Pfam_info" = paste0(Pfam_info, collapse = "; "),
      "Pfam_IPS" = paste0(Pfam_IPS, collapse = "; "),
      "PANTHER_info" = paste0(PANTHER_info, collapse = "; "),
      "PANTHER_IPS" = paste0(PANTHER_IPS, collapse = "; "),
      "multiloc1" = paste0(multiloc1, collapse = "; "),
      "multiloc2" = paste0(multiloc2, collapse = "; "),
      "multiloc3" = paste0(multiloc3, collapse = "; "),
      "KOG_number" = paste0(KOG_number, collapse = "; "),
      "KOG_info" = paste0(KOG_info, collapse = "; "),
      "Keg" = paste0(Keg, collapse = "; "),
      "ko_values" = paste0(ko_values, collapse = "; "),
      "seq_name" = paste0(seq_name, collapse = "; "),
      "UP_AND468" = paste0(UP_AND468, collapse = "; "),
      "UP_OR468" = paste0(UP_OR468, collapse = "; "),
      "DOWN_AND468" = paste0(DOWN_AND468, collapse = "; "),
      "DOWN_OR468" = paste0(DOWN_OR468, collapse = "; "),
      "Increasing" = paste0(Increasing, collapse = "; "),
      "table" = paste0(table, collapse = "; "),
      "set" = paste0(set, collapse = "; "))
} else if (contr == "tx") {
  
  threesets_concat_sig <- rbind(results_anno_onlysig[[1]],
                                results_anno_onlysig[[2]]) |>
    group_by(name_anno) |>
    summarise(
      "RAW.T0-1" =  paste0(`RAW.T0-1`, collapse = "; "),
      "RAW.T0-2" =  paste0(`RAW.T0-2`, collapse = "; "),
      "RAW.T0-3" =  paste0(`RAW.T0-3`, collapse = "; "),
      "RAW.T2-1" =  paste0(`RAW.T2-1`, collapse = "; "),
      "RAW.T2-2" =  paste0(`RAW.T2-2`, collapse = "; "),
      "RAW.T2-3" =  paste0(`RAW.T2-3`, collapse = "; "),
      "RAW.T4-1" =  paste0(`RAW.T4-1`, collapse = "; "),
      "RAW.T4-2" =  paste0(`RAW.T4-2`, collapse = "; "),
      "RAW.T4-3" =  paste0(`RAW.T4-3`, collapse = "; "),
      "RAW.T6-1" =  paste0(`RAW.T6-1`, collapse = "; "),
      "RAW.T6-2" =  paste0(`RAW.T6-2`, collapse = "; "),
      "RAW.T6-3" =  paste0(`RAW.T6-3`, collapse = "; "),
      "RAW.T8-1" =  paste0(`RAW.T8-1`, collapse = "; "),
      "RAW.T8-2" =  paste0(`RAW.T8-2`, collapse = "; "),
      "RAW.T8-3" =  paste0(`RAW.T8-3`, collapse = "; "),
      
      "BATCH_CORRECT.T0-1" =  paste0(`BATCH_CORRECT.T0-1`, collapse = "; "),
      "BATCH_CORRECT.T0-2" =  paste0(`BATCH_CORRECT.T0-2`, collapse = "; "),
      "BATCH_CORRECT.T0-3" =  paste0(`BATCH_CORRECT.T0-3`, collapse = "; "),
      "BATCH_CORRECT.T2-1" =  paste0(`BATCH_CORRECT.T2-1`, collapse = "; "),
      "BATCH_CORRECT.T2-2" =  paste0(`BATCH_CORRECT.T2-2`, collapse = "; "),
      "BATCH_CORRECT.T2-3" =  paste0(`BATCH_CORRECT.T2-3`, collapse = "; "),
      "BATCH_CORRECT.T4-1" =  paste0(`BATCH_CORRECT.T4-1`, collapse = "; "),
      "BATCH_CORRECT.T4-2" =  paste0(`BATCH_CORRECT.T4-2`, collapse = "; "),
      "BATCH_CORRECT.T4-3" =  paste0(`BATCH_CORRECT.T4-3`, collapse = "; "),
      "BATCH_CORRECT.T6-1" =  paste0(`BATCH_CORRECT.T6-1`, collapse = "; "),
      "BATCH_CORRECT.T6-2" =  paste0(`BATCH_CORRECT.T6-2`, collapse = "; "),
      "BATCH_CORRECT.T6-3" =  paste0(`BATCH_CORRECT.T6-3`, collapse = "; "),
      "BATCH_CORRECT.T8-1" =  paste0(`BATCH_CORRECT.T8-1`, collapse = "; "),
      "BATCH_CORRECT.T8-2" =  paste0(`BATCH_CORRECT.T8-2`, collapse = "; "),
      "BATCH_CORRECT.T8-3" =  paste0(`BATCH_CORRECT.T8-3`, collapse = "; "),
      
      "T0_1" = paste0(T0_1, collapse = "; "),
      "T0_2" = paste0(T0_2, collapse = "; "),
      "T0_3" = paste0(T0_3, collapse = "; "),
      "T2_1" = paste0(T2_1, collapse = "; "),
      "T2_2" = paste0(T2_2, collapse = "; "),
      "T2_3" = paste0(T2_3, collapse = "; "),
      "T4_1" = paste0(T4_1, collapse = "; "),
      "T4_2" = paste0(T4_2, collapse = "; "),
      "T4_3" = paste0(T4_3, collapse = "; "),
      "T6_1" = paste0(T6_1, collapse = "; "),
      "T6_2" = paste0(T6_2, collapse = "; "),
      "T6_3" = paste0(T6_3, collapse = "; "),
      "T8_1" = paste0(T8_1, collapse = "; "),
      "T8_2" = paste0(T8_2, collapse = "; "),
      "T8_3" = paste0(T8_3, collapse = "; "),
      
      "T2_vs_T0_ratio" = paste0(T2_vs_T0_ratio, collapse = "; "),
      "T4_vs_T0_ratio" = paste0(T4_vs_T0_ratio, collapse = "; "),
      "T6_vs_T0_ratio" = paste0(T6_vs_T0_ratio, collapse = "; "),
      "T8_vs_T0_ratio" = paste0(T8_vs_T0_ratio, collapse = "; "),
      "T4_vs_T2_ratio" = paste0(T4_vs_T2_ratio, collapse = "; "),
      "T6_vs_T4_ratio" = paste0(T6_vs_T4_ratio, collapse = "; "),
      "T8_vs_T6_ratio" = paste0(T8_vs_T6_ratio, collapse = "; "),
      
      "T2_centered" = paste0(T2_centered, collapse = "; "),
      "T4_centered" = paste0(T4_centered, collapse = "; "),
      "T6_centered" = paste0(T6_centered, collapse = "; "),
      "T8_centered" = paste0(T8_centered, collapse = "; "),
      
      "T2_vs_T0_p.val" = paste0(T2_vs_T0_p.val, collapse = "; "),
      "T4_vs_T0_p.val" = paste0(T4_vs_T0_p.val, collapse = "; "),
      "T6_vs_T0_p.val" = paste0(T6_vs_T0_p.val, collapse = "; "),
      "T8_vs_T0_p.val" = paste0(T8_vs_T0_p.val, collapse = "; "),
      "T4_vs_T2_p.val" = paste0(T4_vs_T2_p.val, collapse = "; "),
      "T6_vs_T4_p.val" = paste0(T6_vs_T4_p.val, collapse = "; "),
      "T8_vs_T6_p.val" = paste0(T8_vs_T6_p.val, collapse = "; "),
      
      "T2_vs_T0_p.adj" = paste0(T2_vs_T0_p.adj, collapse = "; "),
      "T4_vs_T0_p.adj" = paste0(T4_vs_T0_p.adj, collapse = "; "),
      "T6_vs_T0_p.adj" = paste0(T6_vs_T0_p.adj, collapse = "; "),
      "T8_vs_T0_p.adj" = paste0(T8_vs_T0_p.adj, collapse = "; "),
      "T4_vs_T2_p.adj" = paste0(T4_vs_T2_p.adj, collapse = "; "),
      "T6_vs_T4_p.adj" = paste0(T6_vs_T4_p.adj, collapse = "; "),
      "T8_vs_T6_p.adj" = paste0(T8_vs_T6_p.adj, collapse = "; "),
      
      "T2_vs_T0_significant" = paste0(T2_vs_T0_significant, collapse = "; "),
      "T4_vs_T0_significant" = paste0(T4_vs_T0_significant, collapse = "; "),
      "T6_vs_T0_significant" = paste0(T6_vs_T0_significant, collapse = "; "),
      "T8_vs_T0_significant" = paste0(T8_vs_T0_significant, collapse = "; "),
      "T4_vs_T2_significant" = paste0(T4_vs_T2_significant, collapse = "; "),
      "T6_vs_T4_significant" = paste0(T6_vs_T4_significant, collapse = "; "),
      "T8_vs_T6_significant" = paste0(T8_vs_T6_significant, collapse = "; "),
      "significant" = paste0(significant, collapse = "; "),
      
      "T2_vs_T0_CI.L" = paste0(T2_vs_T0_CI.L, collapse = "; "),
      "T2_vs_T0_CI.R" = paste0(T2_vs_T0_CI.R, collapse = "; "),
      "T2_vs_T0_diff" = paste0(T2_vs_T0_diff, collapse = "; "),
      
      "T4_vs_T0_CI.L" = paste0(T4_vs_T0_CI.L, collapse = "; "),
      "T4_vs_T0_CI.R" = paste0(T4_vs_T0_CI.R, collapse = "; "),
      "T4_vs_T0_diff" = paste0(T4_vs_T0_diff, collapse = "; "),
      
      "T6_vs_T0_CI.L" = paste0(T6_vs_T0_CI.L, collapse = "; "),
      "T6_vs_T0_CI.R" = paste0(T6_vs_T0_CI.R, collapse = "; "),
      "T6_vs_T0_diff" = paste0(T6_vs_T0_diff, collapse = "; "),
      
      "T8_vs_T0_CI.L" = paste0(T8_vs_T0_CI.L, collapse = "; "),
      "T8_vs_T0_CI.R" = paste0(T8_vs_T0_CI.R, collapse = "; "),
      "T8_vs_T0_diff" = paste0(T8_vs_T0_diff, collapse = "; "),
      
      "T4_vs_T2_CI.L" = paste0(T4_vs_T2_CI.L, collapse = "; "),
      "T4_vs_T2_CI.R" = paste0(T4_vs_T2_CI.R, collapse = "; "),
      "T4_vs_T2_diff" = paste0(T4_vs_T2_diff, collapse = "; "),
      
      "T6_vs_T4_CI.L" = paste0(T6_vs_T4_CI.L, collapse = "; "),
      "T6_vs_T4_CI.R" = paste0(T6_vs_T4_CI.R, collapse = "; "),
      "T6_vs_T4_diff" = paste0(T6_vs_T4_diff, collapse = "; "),
      
      "T8_vs_T6_CI.L" = paste0(T8_vs_T6_CI.L, collapse = "; "),
      "T8_vs_T6_CI.R" = paste0(T8_vs_T6_CI.R, collapse = "; "),
      "T8_vs_T6_diff" = paste0(T8_vs_T6_diff, collapse = "; "),
      
      "FC_T2_vs_T0" = paste0(FC_T2_vs_T0, collapse = "; "),
      "FC_T4_vs_T0" = paste0(FC_T4_vs_T0, collapse = "; "),
      "FC_T6_vs_T0" = paste0(FC_T6_vs_T0, collapse = "; "),
      "FC_T8_vs_T0" = paste0(FC_T8_vs_T0, collapse = "; "),
      "FC_T4_vs_T2" = paste0(FC_T4_vs_T2, collapse = "; "),
      "FC_T6_vs_T4" = paste0(FC_T6_vs_T4, collapse = "; "),
      "FC_T8_vs_T6" = paste0(FC_T8_vs_T6, collapse = "; "),
      
      "T2v0.reg" = paste0(T2v0.reg, collapse = "; "),
      "T4v0.reg" = paste0(T4v0.reg, collapse = "; "),
      "T6v0.reg" = paste0(T6v0.reg, collapse = "; "),
      "T8v0.reg" = paste0(T8v0.reg, collapse = "; "),
      "T4v2.reg" = paste0(T4v2.reg, collapse = "; "),
      "T6v4.reg" = paste0(T6v4.reg, collapse = "; "),
      "T8v6.reg" = paste0(T8v6.reg, collapse = "; "),
      
      
      "orthocluster" = paste0(orthocluster, collapse = "; "),
      "Aplke1" = paste0(Aplke1, collapse = "; "),
      "Aurli1" = paste0(Aurli1, collapse = "; "),
      "Honfer1" = paste0(Honfer1, collapse = "; "),
      "LabyHa" = paste0(LabyHa, collapse = "; "),
      "Schag1" = paste0(Schag1, collapse = "; "),
      "num.species" = paste0(num.species, collapse = "; "),
      "conserved" = paste0(conserved, collapse = "; "),
      
      "above_median_exp_t0" = paste0(above_median_exp_t0, collapse = "; "),
      "Annotations" = paste0(Annotations, collapse = "; "),
      "genbank" = paste0(genbank, collapse = "; "),
      "Uniprot_Hon" = paste0(Uniprot_Hon, collapse = "; "),
      "gotermId" = paste0(gotermId, collapse = "; "),
      "goName" = paste0(goName, collapse = "; "),
      "gotermType" = paste0(gotermType, collapse = "; "),
      "goAcc" = paste0(goAcc, collapse = "; "),
      "kogdefline" = paste0(kogdefline, collapse = "; "),
      "kogClass" = paste0(kogClass, collapse = "; "),
      "GOcompartment" = paste0(GOcompartment, collapse = "; "),
      "Sig_Pfam_info" = paste0(Sig_Pfam_info, collapse = "; "),
      "Sig_Pfam_number" = paste0(Sig_Pfam_number, collapse = "; "),
      "Pfam_info" = paste0(Pfam_info, collapse = "; "),
      "Pfam_IPS" = paste0(Pfam_IPS, collapse = "; "),
      "PANTHER_info" = paste0(PANTHER_info, collapse = "; "),
      "PANTHER_IPS" = paste0(PANTHER_IPS, collapse = "; "),
      "multiloc1" = paste0(multiloc1, collapse = "; "),
      "multiloc2" = paste0(multiloc2, collapse = "; "),
      "multiloc3" = paste0(multiloc3, collapse = "; "),
      "KOG_number" = paste0(KOG_number, collapse = "; "),
      "KOG_info" = paste0(KOG_info, collapse = "; "),
      "Keg" = paste0(Keg, collapse = "; "),
      "ko_values" = paste0(ko_values, collapse = "; "),
      "seq_name" = paste0(seq_name, collapse = "; "),
      "UP_AND468" = paste0(UP_AND468, collapse = "; "),
      "UP_OR468" = paste0(UP_OR468, collapse = "; "),
      "DOWN_AND468" = paste0(DOWN_AND468, collapse = "; "),
      "DOWN_OR468" = paste0(DOWN_OR468, collapse = "; "),
      "UP_T6T4_OR_T4T2" = paste0(UP_T6T4_OR_T4T2, collapse = "; "),
      "Increasing" = paste0(Increasing, collapse = "; "),
      "table" = paste0(table, collapse = "; "),
      "set" = paste0(set, collapse = "; "))
} else {
  # Code to execute if 'set' is neither "t0" nor "tx"
  # Add your code here
}

threesets_concat |> 
  write_csv(file.path("proteomics/DEP_results_IndependentSets", paste0("df_anno3s.",contr,".csv")))

threesets_concat_sig |> 
  write_csv(file.path("proteomics/DEP_results_IndependentSets", paste0("df_anno3s_filtsig.",contr,".csv")))

## Export significant proteins as fasta files

# Fasta of all significant
cand1_all <- results_anno |> map(ungroup)
cand1_sig <- cand1_all |>
  map(filter, significant == TRUE) |>
  map(select,seq_name,sequence)

for (j in seq_along(cand1_all)) {
  d <- do.call(rbind, lapply(seq(nrow(cand1_sig[[j]])),
                             function(i) t(cand1_sig[[j]][i, ])))
  write.table(d,
              file.path("proteomics/DEP_results_fasta_IndependentSets", paste0("cand_sig.", file_set[j],".", contr,".fasta")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Fasta of all up
cand1_all <- results_anno |> map(ungroup)
cand1_up <- cand1_all |>
  map(filter, T2v0.reg == "UP" | T4v0.reg == "UP" |
        T6v0.reg == "UP" | T8v0.reg == "UP") |>
  map(select,seq_name,sequence)

for (j in seq_along(cand1_up)) {
  d <- do.call(rbind, lapply(seq(nrow(cand1_up[[j]])),
                             function(i) t(cand1_up[[j]][i, ])))
  write.table(d,
              file.path("proteomics/DEP_results_fasta_IndependentSets", paste0("cand_up.", file_set[j],".", contr,".fasta")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Fasta of all down
cand1_all <- results_anno |> map(ungroup)
cand1_down <- cand1_all |>
  map(filter, T2v0.reg == "DOWN" | T4v0.reg == "DOWN" |
        T6v0.reg == "DOWN" | T8v0.reg == "DOWN") |>
  map(select,seq_name,sequence)


for (j in seq_along(cand1_down)) {
  d <- do.call(rbind, lapply(seq(nrow(cand1_down[[j]])),
                             function(i) t(cand1_down[[j]][i, ])))
  write.table(d,
              file.path("proteomics/DEP_results_fasta_IndependentSets", paste0("cand_down.", file_set[j],".", contr,".fasta")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}


## Volcano plots colored by KOGclass

# Volcano plot for annotated set.
results_anno[[3]] |> colnames()

df <- results_anno[[1]] |> 
  select(name_anno,ID,significant,cluster,kogClass,ko_values,
         T2_vs_T0_ratio,T4_vs_T0_ratio,T6_vs_T0_ratio,T8_vs_T0_ratio,
         T2_vs_T0_p.val,T4_vs_T0_p.val,T6_vs_T0_p.val,T8_vs_T0_p.val) |> 
  rename(Sig.Any = significant) |> 
  mutate(kogClass=ifelse(is.na(kogClass),"Unknown kogClass",kogClass)) |> 
  separate(kogClass,sep="; ",into=LETTERS[1:4]) |> 
  pivot_longer(names_to = "kogClass_name",
               values_to = "kogClass",5:8) |>    #kogClass
  filter(!is.na(kogClass)) |> 
  mutate(ko_values=ifelse(is.na(ko_values),"Unknown ko_value",ko_values)) |> 
  separate(ko_values,sep="; ",into=LETTERS[1:20]) |> 
  pivot_longer(names_to = "ko_values_name",
               values_to = "ko_values",5:24) |> #ko_values
  filter(!is.na(ko_values)) |> 
  select(-c(ko_values_name,kogClass_name)) |>  
  pivot_longer(names_to = "type", values_to = "values1", 
               c(5:8, #ratios
                 9:12)) |>  #p_values
  separate(type, sep = "_", into = c("C1", "C2", "C3", "type")) |> 
  mutate(Time = paste0(C1, "_", C2, "_", C3)) |> 
  select(-c("C1", "C2", "C3", "ID")) |> 
  pivot_wider(names_from = type, values_from = values1)|> 
  filter(Time %in% c("T2_vs_T0","T4_vs_T0","T6_vs_T0","T8_vs_T0")) |> 
  unchop(everything()) |> 
  mutate(neg.log.p.val = -log(p.val)) |>
  mutate(Time = ifelse(Time == "T2_vs_T0", "T2 vs T0",
                       ifelse(Time == "T4_vs_T0", "T4 vs T0",
                              ifelse(Time == "T6_vs_T0", "T6 vs T0",
                                     ifelse(Time == "T8_vs_T0", "T8 vs T0", "NA")))))

## KogClass colored by clusters.
list_plot <- list()
for (i in 1:length(unique(df$kogClass))){
  vec <- unique(df$kogClass) |> sort()
  kog <- vec[i]
  print(kog)
  
  data_kog <- df|>
    filter(kogClass==kog)
  
  data_sig_false <- data_kog|>
    filter(cluster=="not significant")
  data_sig_true <- data_kog|>
    filter(cluster!="not significant")
  
  thr <- df |> 
    filter(Sig.Any==FALSE) |>
    group_by(Time) |>
    summarise(max.nlp = max(neg.log.p.val))
  
  list_plot[[i]] <- ggplot() +
    geom_point(aes(x = ratio, y = neg.log.p.val,
                   col = cluster),
               data = data_sig_false)+
    geom_point(aes(x = ratio, y = neg.log.p.val,
                   col = cluster),
               data = data_sig_true) +
    facet_wrap(.~Time) +
    geom_hline(data = thr, aes(yintercept = max.nlp),linetype='dotted') +
    geom_vline(xintercept=0)+
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_color_manual(breaks = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4","not significant"),values = c(brewer.pal(5, "Set1"), "grey")) +
    xlab(expression('log'[2]*'(Fold change)'))+
    ylab(expression('-log'[10]*'(p value)'))+
    ggtitle(kog)+
    xlim(-1.5, +1.5)+
    ylim(-5,35)+
    geom_text_repel(aes(ratio, neg.log.p.val, label = name_anno),size = 1,data=data_sig_true)+
    geom_vline(xintercept=log2(fold_threshold),linetype='dotted')+
    geom_vline(xintercept=log2(inverse_fold),linetype='dotted')
}

## KogClass colored by clusters.
pdf(paste0("proteomics/img_IndependentSets/KogClass_percluster.",contr,".pdf"))
list_plot
dev.off()


## KogClass colored by ko_values
list_plot <- list()
for (i in 25){
  vec <- unique(df$kogClass) |> sort()
  kog <- vec[i]
  print(kog)
  
  data_kog <- df|>
    filter(kogClass==kog)
  
  data_sig_false <- data_kog|>
    filter(!ko_values%in%c('ko03011',"ko03009","ko03016","ko03012"))|> 
    mutate(ko_values=ifelse(ko_values=='ko03011',"Ribosome",
                            ifelse(ko_values=='ko03009',"Ribosome biogenesis",
                                   ifelse(ko_values=="ko03016","Transfer RNA biogenesis",
                                          ifelse(ko_values=="ko03012","Translation factors",
                                                 ifelse(ko_values=="other","other",ko_values))))))
  
  data_sig_true <- data_kog|>
    filter(ko_values%in%c('ko03011',"ko03009","ko03016","ko03012")) |> 
    mutate(ko_values=ifelse(ko_values=='ko03011',"Ribosome",
                            ifelse(ko_values=='ko03009',"Ribosome biogenesis",
                                   ifelse(ko_values=="ko03016","Transfer RNA biogenesis",
                                          ifelse(ko_values=="ko03012","Translation factors",
                                                 ifelse(ko_values=="other","other",ko_values))))))
  
  thr <- df |> 
    filter(Sig.Any==FALSE) |>
    group_by(Time) |>
    summarise(max.nlp = max(neg.log.p.val))
  
  list_plot[[i]] <- ggplot() +
    geom_point(aes(x = ratio, y = neg.log.p.val,
                   col = "other"),
               data = data_sig_false)+
    geom_point(aes(x = ratio, y = neg.log.p.val,
                   col = ko_values),
               data = data_sig_true) +
    facet_wrap(.~Time) +
    geom_hline(data = thr, aes(yintercept = max.nlp),linetype='dotted') +
    geom_vline(xintercept=0)+
    theme_minimal() +
    theme(legend.position = "bottom") +
    xlab(expression('log'[2]*'(Fold change)'))+
    ylab(expression('-log'[10]*'(p value)'))+
    ggtitle(kog)+
    xlim(-1.10, +1)+
    ylim(-5,35)+
    scale_color_manual(breaks = c("Ribosome biogenesis", "Ribosome", "Translation factors", "Transfer RNA biogenesis","other"),values = c(brewer.pal(4, "Set1"), "grey"))+
    geom_hline(yintercept = -log(0.05,10))+
    geom_vline(xintercept=log2(fold_threshold),linetype='dotted')+
    geom_vline(xintercept=log2(inverse_fold),linetype='dotted')
  
  
  # geom_text_repel(aes(ratio, neg.log.p.val, label = name_anno),size = 1,data=data_sig_true)
}

## KogClass colored by ko_values
pdf(paste0("proteomics/img_IndependentSets/Transcription_plot.",contr,".pdf"))
list_plot[[25]]
dev.off()

png_path <- paste0("proteomics/img_IndependentSets/Transcription_plot.",contr,".png")
png(png_path,
    width = 4 * a,
    height = 4 * a,
    units = "in",
    res = 300)
list_plot[[25]]
dev.off()

## Fisher's tests. Only do this for contr='t0'

## get_table_test(class) class:kogClass/ko_values

#get_table_test("kogClass")

get_table_test <- function(class){
  
# ALLDET
table_m <- results_anno[[3]] |> 
    separate_rows(class, sep = "; ")  |> 
    mutate(cat=ifelse(is.na(eval(parse(text = class))),"NA",eval(parse(text = class)))) |>    
    filter(cat!="NA") |> # This removes the proteins that are NA
    select(name_anno,significant,cluster,cat)
  
#Count all significant in ALLDET
table_s <- results_anno[[3]] |> 
    separate_rows(class, sep = "; ")  |> 
    mutate(cat=ifelse(is.na(eval(parse(text = class))),"NA",eval(parse(text = class)))) |>    
    filter(cat!="NA") |> # This removes the proteins that are NA
    select(name_anno,significant,cluster,cat)
  
# Count all detected in all sets.
table_all <- rbind(results_anno[[1]],results_anno[[2]]) |> 
    select(name_anno,significant,cluster,class)|> 
    separate_rows(class, sep = "; ")  |> 
    mutate(cat=ifelse(is.na(eval(parse(text = class))),"NA",eval(parse(text = class)))) |>    
    filter(cat!="NA") |> 
    select(name_anno,significant,cluster,cat)
  
# Unique proteins detected across all non NA categories.
total_detected_m <- table_m |> select(name_anno) |> pull() |>  unique() |> length()
total_sig_m <- table_m |> filter(significant==TRUE) |> select(  name_anno) |> pull() |> unique() |> length()
  total_down_m <- table_m |> filter(significant==TRUE) |> 
    filter(cluster=="Cluster 1"|cluster=="Cluster 2")|> select(name_anno) |>  pull() |>unique() |> length()
  total_up_m <- table_m |> filter(significant==TRUE) |> 
    filter(cluster=="Cluster 3"|cluster=="Cluster 4")|> select(name_anno) |>pull() |> unique() |> length()
  total_clu1_m <- table_m |> filter(significant==TRUE) |> 
    filter(cluster=="Cluster 1")|> select(name_anno) |> pull() |>unique() |> length()
  total_clu2_m <- table_m |> filter(significant==TRUE) |> 
    filter(cluster=="Cluster 2")|> select(name_anno) |> pull() |>unique() |> length()
  total_clu3_m <- table_m |> filter(significant==TRUE) |> 
    filter(cluster=="Cluster 3")|> select(name_anno) |> pull() |>unique() |> length()
  total_clu4_m <- table_m |> filter(significant==TRUE) |> 
    filter(cluster=="Cluster 4")|> select(name_anno) |> pull() |>unique() |> length()

  total_detected_s <- table_all |> select(name_anno) |> pull() |>unique() |> length()
  total_sig_s <- table_s |> filter(significant==TRUE) |> select(  name_anno) |> pull() |> unique() |> length()
  total_down_s <- table_s |> filter(significant==TRUE) |> 
    filter(cluster=="Cluster 1"|cluster=="Cluster 2")|> select(name_anno) |>  pull() |>unique() |> length()
  total_up_s <- table_s |> filter(significant==TRUE) |> 
    filter(cluster=="Cluster 3"|cluster=="Cluster 4")|> select(name_anno) |>pull() |> unique() |> length()
  total_clu1_s <- table_s |> filter(significant==TRUE) |> 
    filter(cluster=="Cluster 1")|> select(name_anno) |> pull() |>unique() |> length()
  total_clu2_s <- table_s |> filter(significant==TRUE) |> 
    filter(cluster=="Cluster 2")|> select(name_anno) |> pull() |>unique() |> length()
  total_clu3_s <- table_s |> filter(significant==TRUE) |> 
    filter(cluster=="Cluster 3")|> select(name_anno) |> pull() |>unique() |> length()
  total_clu4_s <- table_s |> filter(significant==TRUE) |> 
    filter(cluster=="Cluster 4")|> select(name_anno) |> pull() |>unique() |> length()

  cat_vect_m <-table_m |> 
    select(cat)|> distinct() |> pull() |> sort()
  cat_vect_s <-table_s |> 
    select(cat)|> distinct() |> pull() |> sort()
  cat_vect <- c(cat_vect_m,cat_vect_s) |> unique()
  
  df_fisher <- list()
  for (i in cat_vect){
    # For every category:
    
    # Count how many unique proteins were detected in that category.
    detected_val_m <- table_m |> 
      filter(cat==i) |> 
      select(name_anno) |> pull() |> unique() |> length()
    detected_val_s <- table_all |> 
      filter(cat==i) |> 
      select(name_anno) |> pull() |> unique() |> length()
    
    # Count how many unique proteins were significant in that category.
    significant_val_m <- table_m |> 
      filter(cat==i) |> 
      filter(significant==TRUE) |> 
      select(name_anno) |> pull() |> unique() |> length()
    significant_val_s <- table_s |> 
      filter(cat==i) |> 
      filter(significant==TRUE) |> 
      select(name_anno) |> pull() |> unique() |> length()
    
    
    # Count how many unique proteins were in cl_X in that category.
    clu1_val_m <- table_m |> 
      filter(cat==i) |> 
      filter(significant==TRUE) |> 
      filter(cluster=="Cluster 1") |> 
      select(name_anno) |> pull() |> unique() |> length()
    clu1_val_s <- table_s |> 
      filter(cat==i) |> 
      filter(significant==TRUE) |> 
      filter(cluster=="Cluster 1") |> 
      select(name_anno) |> pull() |> unique() |> length()
    clu2_val_m <- table_m |> 
      filter(cat==i) |> 
      filter(significant==TRUE) |> 
      filter(cluster=="Cluster 2") |> 
      select(name_anno) |> pull() |> unique() |> length()
    clu2_val_s <- table_s |> 
      filter(cat==i) |> 
      filter(significant==TRUE) |> 
      filter(cluster=="Cluster 2") |> 
      select(name_anno) |> pull() |> unique() |> length()
    clu3_val_m <- table_m |> 
      filter(cat==i) |> 
      filter(significant==TRUE) |> 
      filter(cluster=="Cluster 3") |> 
      select(name_anno) |> pull() |> unique() |> length()
    clu3_val_s <- table_s |> 
      filter(cat==i) |> 
      filter(significant==TRUE) |> 
      filter(cluster=="Cluster 3") |> 
      select(name_anno) |> pull() |> unique() |> length()
    clu4_val_m <- table_m |> 
      filter(cat==i) |> 
      filter(significant==TRUE) |> 
      filter(cluster=="Cluster 4") |> 
      select(name_anno) |> pull() |> unique() |> length()
    clu4_val_s <- table_s |> 
      filter(cat==i) |> 
      filter(significant==TRUE) |> 
      filter(cluster=="Cluster 4") |> 
      select(name_anno) |> pull() |> unique() |> length()
    
    # How many unique proteins were in clusters 1 or 2 in that category
    down_m <- table_m |> 
      filter(cat==i) |> 
      filter(significant==TRUE) |> 
      filter(cluster=="Cluster 1"|cluster=="Cluster 2") |> 
      select(name_anno) |> pull() |> unique() |> length()
    down_s <- table_s |> 
      filter(cat==i) |> 
      filter(significant==TRUE) |> 
      filter(cluster=="Cluster 1"|cluster=="Cluster 2") |> 
      select(name_anno) |> pull() |> unique() |> length()
    
    up_m <- table_m |> 
      filter(cat==i) |> 
      filter(significant==TRUE) |> 
      filter(cluster=="Cluster 3"|cluster=="Cluster 4") |> 
      select(name_anno) |> pull() |> unique() |> length()
    up_s <- table_s |> 
      filter(cat==i) |> 
      filter(significant==TRUE) |> 
      filter(cluster=="Cluster 3"|cluster=="Cluster 4") |> 
      select(name_anno) |> pull() |> unique() |> length()
    
    df_fisher[[i]] <- c(detected_val_m,significant_val_m,down_m,up_m,clu1_val_m,clu2_val_m,clu3_val_m,clu4_val_m,
                        detected_val_s,significant_val_s,down_s,up_s,clu1_val_s,clu2_val_s,clu3_val_s,clu4_val_s)
  }
  df_fisher_df <- do.call(rbind, df_fisher)
  df_fisher_df <- cbind(cat_vect,df_fisher_df)
  
  rownames(df_fisher_df) <- NULL
  colnames(df_fisher_df) <- c("category","detected_m","significant_m",
                              "down_m","up_m","clu1_m","clu2_m","clu3_m",
                              "clu4_m","detected_s","significant_s","down_s",
                              "up_s","clu1_s","clu2_s","clu3_s","clu4_s")
  
  df_fisher_df <- df_fisher_df |>
    as.data.frame()
  
  df_fisher_df$total_detected_m <- total_detected_m
  df_fisher_df$total_sig_m <- total_sig_m
  df_fisher_df$total_down_m <- total_down_m
  df_fisher_df$total_up_m <- total_up_m
  df_fisher_df$total_clu1_m <- total_clu1_m
  df_fisher_df$total_clu2_m <- total_clu2_m
  df_fisher_df$total_clu3_m <- total_clu3_m
  df_fisher_df$total_clu4_m <- total_clu4_m

  df_fisher_df$total_detected_s <- total_detected_s
  df_fisher_df$total_sig_s <- total_sig_s
  df_fisher_df$total_down_s <- total_down_s
  df_fisher_df$total_up_s <- total_up_s
  df_fisher_df$total_clu1_s <- total_clu1_s
  df_fisher_df$total_clu2_s <- total_clu2_s
  df_fisher_df$total_clu3_s <- total_clu3_s
  df_fisher_df$total_clu4_s <- total_clu4_s

  df_fisher_df <- mutate_at(df_fisher_df, 2:33, as.numeric)
  
  sig_m <- c()
  down_m <- c()
  up_m <- c()
  clu1_m <- c()
  clu2_m <- c()
  clu3_m <- c()
  clu4_m <- c()

  sig_s <- c()
  down_s <- c()
  up_s <- c()
  clu1_s <- c()
  clu2_s <- c()
  clu3_s <- c()
  clu4_s <- c()

  rows <- dim(df_fisher_df)[1]
  for (i in 1:rows){
    b <- df_fisher_df[i,]
    
    sig_m[i] <- dhyper(b$significant_m,b$total_sig_m, b$total_detected_m-b$total_sig_m,b$detected_m)
    down_m[i] <- dhyper(b$down_m,b$total_down_m, b$total_sig_m-b$total_down_m,b$significant_m)
    up_m[i] <- dhyper(b$up_m,b$total_up_m, b$total_sig_m-b$total_up_m,b$significant_m)
    clu1_m[i] <- dhyper(b$clu1_m,b$total_clu1_m, b$total_sig_m-b$total_clu1_m,b$significant_m)
    clu2_m[i] <- dhyper(b$clu2_m,b$total_clu2_m, b$total_sig_m-b$total_clu2_m,b$significant_m)
    clu3_m[i] <- dhyper(b$clu3_m,b$total_clu3_m, b$total_sig_m-b$total_clu3_m,b$significant_m)
    clu4_m[i] <- dhyper(b$clu4_m,b$total_clu4_m, b$total_sig_m-b$total_clu4_m,b$significant_m)
    
    sig_s[i] <- dhyper(b$significant_s,b$total_sig_s, b$total_detected_s-b$total_sig_s,b$detected_s)
    down_s[i] <- dhyper(b$down_s,b$total_down_s, b$total_sig_s-b$total_down_s,b$significant_s)
    up_s[i] <- dhyper(b$up_s,b$total_up_s, b$total_sig_s-b$total_up_s,b$significant_s)
    clu1_s[i] <- dhyper(b$clu1_s,b$total_clu1_s, b$total_sig_s-b$total_clu1_s,b$significant_s)
    clu2_s[i] <- dhyper(b$clu2_s,b$total_clu2_s, b$total_sig_s-b$total_clu2_s,b$significant_s)
    clu3_s[i] <- dhyper(b$clu3_s,b$total_clu3_s, b$total_sig_s-b$total_clu3_s,b$significant_s)
    clu4_s[i] <- dhyper(b$clu4_s,b$total_clu4_s, b$total_sig_s-b$total_clu4_s,b$significant_s)

    
  }
  
  df_fisher_df$p.val_sig_vs_det_m <- sig_m
  df_fisher_df$p.val_down_vs_sig_m <- down_m
  df_fisher_df$p.val_up_vs_sig_m <- up_m
  df_fisher_df$p.val_clu1_vs_sig_m <- clu1_m
  df_fisher_df$p.val_clu2_vs_sig_m <- clu2_m
  df_fisher_df$p.val_clu3_vs_sig_m <- clu3_m
  df_fisher_df$p.val_clu4_vs_sig_m <- clu4_m
  df_fisher_df$p.val_sig_vs_det_m <- sig_m
  
  df_fisher_df$p.val_sig_vs_det_s <- sig_s
  df_fisher_df$p.val_down_vs_sig_s <- down_s
  df_fisher_df$p.val_up_vs_sig_s <- up_s
  df_fisher_df$p.val_clu1_vs_sig_s <- clu1_s
  df_fisher_df$p.val_clu2_vs_sig_s <- clu2_s
  df_fisher_df$p.val_clu3_vs_sig_s <- clu3_s
  df_fisher_df$p.val_clu4_vs_sig_s <- clu4_s

  df_fisher_df_res <- df_fisher_df |> 
    mutate(p.adj_sig_vs_det_m=p.adjust(p.val_sig_vs_det_m,
                                       method = "BH"),
           p.adj_down_vs_sig_m=p.adjust(p.val_down_vs_sig_m,
                                        method = "BH"),
           p.adj_up_vs_sig_m=p.adjust(p.val_up_vs_sig_m,
                                      method = "BH"),
           p.adj_clu1_vs_sig_m=p.adjust(p.val_clu1_vs_sig_m,
                                        method = "BH"),
           p.adj_clu2_vs_sig_m=p.adjust(p.val_clu2_vs_sig_m,
                                        method = "BH"),
           p.adj_clu3_vs_sig_m=p.adjust(p.val_clu3_vs_sig_m,
                                        method = "BH"),
           p.adj_clu4_vs_sig_m=p.adjust(p.val_clu4_vs_sig_m,
                                        method = "BH")) |> 
    mutate(p.adj_sig_vs_det_s=p.adjust(p.val_sig_vs_det_s,
                                       method = "BH"),
           p.adj_down_vs_sig_s=p.adjust(p.val_down_vs_sig_s,
                                        method = "BH"),
           p.adj_up_vs_sig_s=p.adjust(p.val_up_vs_sig_s,
                                      method = "BH"),
           p.adj_clu1_vs_sig_s=p.adjust(p.val_clu1_vs_sig_s,
                                        method = "BH"),
           p.adj_clu2_vs_sig_s=p.adjust(p.val_clu2_vs_sig_s,
                                        method = "BH"),
           p.adj_clu3_vs_sig_s=p.adjust(p.val_clu3_vs_sig_s,
                                        method = "BH"),
           p.adj_clu4_vs_sig_s=p.adjust(p.val_clu4_vs_sig_s,
                                        method = "BH"))
  
  res <- df_fisher_df_res |> mutate(test_sig_vs_det_m=ifelse(p.adj_sig_vs_det_m>0.05,
                                                             "NS",ifelse(significant_m/detected_m>total_sig_m/total_detected_m,"OVER","UNDER")))|> 
    mutate(test_down_vs_sig_m=ifelse(p.adj_down_vs_sig_m>0.05,
                                     "NS",ifelse(down_m/significant_m>total_down_m/total_sig_m,"OVER","UNDER"))) |>  
    mutate(test_up_vs_sig_m=ifelse(p.adj_up_vs_sig_m>0.05,
                                   "NS",ifelse(up_m/significant_m>total_up_m/total_sig_m,"OVER","UNDER"))) |>
    
    mutate(test_clu1_vs_sig_m=ifelse(p.adj_clu1_vs_sig_m>0.05,
                                     "NS",ifelse(clu1_m/significant_m>total_clu1_m/total_sig_m,"OVER","UNDER")))|> 
    mutate(test_clu2_vs_sig_m=ifelse(p.adj_clu2_vs_sig_m>0.05,
                                     "NS",ifelse(clu2_m/significant_m>total_clu2_m/total_sig_m,"OVER","UNDER")))|> 
    mutate(test_clu3_vs_sig_m=ifelse(p.adj_clu3_vs_sig_m>0.05,
                                     "NS",ifelse(clu3_m/significant_m>total_clu3_m/total_sig_m,"OVER","UNDER")))|> 
    mutate(test_clu4_vs_sig_m=ifelse(p.adj_clu4_vs_sig_m>0.05,
                                     "NS",ifelse(clu4_m/significant_m>total_clu4_m/total_sig_m,"OVER","UNDER")))|> 
    mutate(FIGURE_m=ifelse(test_sig_vs_det_m=="NS"&test_down_vs_sig_m=="NS"&test_up_vs_sig_m=="NS","Not in figure","In figure")) |> 
    mutate(test_sig_vs_det_s=ifelse(p.adj_sig_vs_det_s>0.05,
                                    "NS",ifelse(significant_s/detected_s>total_sig_s/total_detected_s,"OVER","UNDER")))|> 
    mutate(test_down_vs_sig_s=ifelse(p.adj_down_vs_sig_s>0.05,
                                     "NS",ifelse(down_s/significant_s>total_down_s/total_sig_s,"OVER","UNDER"))) |>  
    mutate(test_up_vs_sig_s=ifelse(p.adj_up_vs_sig_s>0.05,
                                   "NS",ifelse(up_s/significant_s>total_up_s/total_sig_s,"OVER","UNDER"))) |>
    
    mutate(test_clu1_vs_sig_s=ifelse(p.adj_clu1_vs_sig_s>0.05,
                                     "NS",ifelse(clu1_s/significant_s>total_clu1_s/total_sig_s,"OVER","UNDER")))|> 
    mutate(test_clu2_vs_sig_s=ifelse(p.adj_clu2_vs_sig_s>0.05,
                                     "NS",ifelse(clu2_s/significant_s>total_clu2_s/total_sig_s,"OVER","UNDER")))|> 
    mutate(test_clu3_vs_sig_s=ifelse(p.adj_clu3_vs_sig_s>0.05,
                                     "NS",ifelse(clu3_s/significant_s>total_clu3_s/total_sig_s,"OVER","UNDER")))|> 
    mutate(test_clu4_vs_sig_s=ifelse(p.adj_clu4_vs_sig_s>0.05,
                                     "NS",ifelse(clu4_s/significant_s>total_clu4_s/total_sig_s,"OVER","UNDER")))|> 
    mutate(FIGURE_s=ifelse(test_sig_vs_det_s=="NS"&test_down_vs_sig_s=="NS"&test_up_vs_sig_s=="NS","Not in figure","In figure")) 
  
  res
}

d_kC <- get_table_test("kogClass") |> 
  select(category,detected_m,significant_m,clu1_m,clu2_m,clu3_m,clu4_m,down_m,up_m,total_detected_m,total_sig_m,total_down_m,total_up_m,total_clu1_m,total_clu2_m,total_clu3_m,total_clu4_m,p.val_sig_vs_det_m,p.val_down_vs_sig_m,p.val_up_vs_sig_m,p.val_clu1_vs_sig_m,p.val_clu2_vs_sig_m,p.val_clu3_vs_sig_m,p.val_clu4_vs_sig_m,p.adj_sig_vs_det_m,p.adj_down_vs_sig_m,p.adj_up_vs_sig_m,p.adj_clu1_vs_sig_m,p.adj_clu2_vs_sig_m,p.adj_clu3_vs_sig_m,p.adj_clu4_vs_sig_m,test_sig_vs_det_m,test_down_vs_sig_m,test_up_vs_sig_m,test_clu1_vs_sig_m,test_clu2_vs_sig_m,test_clu3_vs_sig_m,test_clu4_vs_sig_m,FIGURE_m)
d_kv <- get_table_test("ko_values")|> 
  select(category,detected_m,significant_m,clu1_m,clu2_m,clu3_m,clu4_m,down_m,up_m,total_detected_m,total_sig_m,total_down_m,total_up_m,total_clu1_m,total_clu2_m,total_clu3_m,total_clu4_m,p.val_sig_vs_det_m,p.val_down_vs_sig_m,p.val_up_vs_sig_m,p.val_clu1_vs_sig_m,p.val_clu2_vs_sig_m,p.val_clu3_vs_sig_m,p.val_clu4_vs_sig_m,p.adj_sig_vs_det_m,p.adj_down_vs_sig_m,p.adj_up_vs_sig_m,p.adj_clu1_vs_sig_m,p.adj_clu2_vs_sig_m,p.adj_clu3_vs_sig_m,p.adj_clu4_vs_sig_m,test_sig_vs_det_m,test_down_vs_sig_m,test_up_vs_sig_m,test_clu1_vs_sig_m,test_clu2_vs_sig_m,test_clu3_vs_sig_m,test_clu4_vs_sig_m,FIGURE_m)

write_csv(d_kC,paste0("proteomics/DEP_results_IndependentSets/","ET_data_kogClass",".csv"))
write_csv(d_kv,paste0("proteomics/DEP_results_IndependentSets/","ET_data_ko_values",".csv"))


# Reorganize manually to make input for figure:

## Fisher's plot

FT <- list()
#FT[[1]] <- read_csv("proteomics/DEP_results_MultiMapping/ET_data_FisherTables_Nov2023.csv")
#FT[[2]] <- read_csv("proteomics/DEP_results_IndependentSets/ET_data_FisherTables_March2023_Fisher_fig_s.csv")
FT[[1]] <- read_csv("proteomics/DEP_results_IndependentSets/ET_data_FisherTables_March2023 - Fisher_fig_IndSets.csv")

colnames(FT[[1]]) <- c("category","detected","significant",
                       "clu1","clu2","clu3","clu4","down","up",
                       "test_sig_vs_det","test_up_vs_sig",
                       "test_clu1_vs_sig","test_clu2_vs_sig","test_clu3_vs_sig","test_clu4_vs_sig",
                       "p.adj_sig_vs_det","p.adj_up_vs_sig",
                       "p.adj_clu1_vs_sig","p.adj_clu2_vs_sig","p.adj_clu3_vs_sig","p.adj_clu4_vs_sig")

#colnames(FT[[2]]) <- c("category","detected","significant","clu1","clu2","clu3","clu4","down","up","test_sig_vs_det","test_up_vs_sig","test_clu1_vs_sig","test_clu2_vs_sig","test_clu3_vs_sig","test_clu4_vs_sig","p.adj_sig_vs_det","p.adj_up_vs_sig","p.adj_clu1_vs_sig","p.adj_clu2_vs_sig","p.adj_clu3_vs_sig","p.adj_clu4_vs_sig")

numb <- 1

FT1 <- FT[[1]] |> 
  mutate(category=factor(category,
                         levels=FT[[1]]$category))


# New Fisher table
panel1 <- FT1  |> 
  ggplot() +
  geom_bar(aes(x = -1*down, 
               y = category,fill="#f1a340"),stat="identity")+
  geom_bar(aes(x = up, 
               y = category,fill="#998ec3ff"),stat="identity")+
  theme_minimal()+
  theme(legend.position = "none") +
  ylab("") +
  xlab("Number of significant proteins") +
  geom_vline(xintercept = 0) +
  scale_y_discrete(limits = rev)+
  geom_text(aes(x=50,y=category,label=detected))+
  scale_x_continuous(breaks=c(-30,-20,-10,0,10,20,30,40))+
  xlim(-35,50)+
  theme(text = element_text(size = 20))  

small_FT <- FT1 |> 
  select(c(1,starts_with("test"))) |> 
  mutate(SigDE=as.factor(test_sig_vs_det),
         Up_Down=as.factor(test_up_vs_sig),
         CLU1_vs_SIG=as.factor(test_clu1_vs_sig),
         CLU2_vs_SIG=as.factor(test_clu2_vs_sig),
         CLU3_vs_SIG=as.factor(test_clu3_vs_sig),
         CLU4_vs_SIG=as.factor(test_clu4_vs_sig)) |> 
  pivot_longer(names_to="test",values_to = "result",2:7)  |>
  mutate(result=factor(result,levels=c("UNDER","NS","OVER")))

#Under 24, NS 20, Over 25 
b <- ggplot()+
  geom_point(aes(x=test,
                 y=category,
                 shape=result,size=result,fill=result),data=small_FT)+
  scale_shape_manual(values=c(25,20,24))+
  scale_size_manual(values=c(6,0,6))+
  scale_fill_manual(values=c('#FFFFFF','#808080','#000000'))+
  theme_minimal()+
  theme(legend.position="right",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_discrete(limits = rev)+
  theme(text = element_text(size = 20)) 

fig7 <- wrap_plots(panel1,b,widths = c(1,0.2))

ggsave(filename = "proteomics/img_IndependentSets/indset.fig7.s.pdf",
       plot=fig7,
       width=20,
       height=10)

ggsave(filename = "proteomics/img_IndependentSets/indset.fig7.s.png",
       plot=fig7,
       width=20,
       height=10,dpi=600)

# Only Sig Det and Up Sig

small_FT <- FT1 |> 
  select(c(1,starts_with("test"))) |> 
  mutate(SigDE=as.factor(test_sig_vs_det),
         Up_Down=as.factor(test_up_vs_sig),
         CLU1_vs_SIG=as.factor(test_clu1_vs_sig),
         CLU2_vs_SIG=as.factor(test_clu2_vs_sig),
         CLU3_vs_SIG=as.factor(test_clu3_vs_sig),
         CLU4_vs_SIG=as.factor(test_clu4_vs_sig)) |> 
  pivot_longer(names_to="test",values_to = "result",2:7)  |>
  mutate(result=factor(result,levels=c("UNDER","NS","OVER"))) |> 
  select(category,test,result) |> 
  filter(!grepl("clu", test))

c <- ggplot()+
  geom_point(aes(x=test,
                 y=category,
                 shape=result,size=result,fill=result),data=small_FT)+
  scale_shape_manual(values=c(25,20,24))+
  scale_size_manual(values=c(6,0,6))+
  scale_fill_manual(values=c('#FFFFFF','#808080','#000000'))+
  theme_minimal()+
  theme(legend.position="right",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_discrete(limits = rev)+
  theme(text = element_text(size = 20)) 

fig7.2 <- wrap_plots(panel1,c,widths = c(1,0.2))

ggsave(filename = "proteomics/img_IndependentSets/indset.fig7.pdf",
       plot=fig7.2,
       width=20,
       height=10)

ggsave(filename = "proteomics/img_IndependentSets/indset.fig7.png",
       plot=fig7.2,
       width=20,
       height=10,dpi=600)


# Check Venn
c12 <- results_anno[[1]] |> filter(cluster=="Cluster 1"|
                              cluster=="Cluster 2") |> 
  pull(Protein_id)

c34 <- results_anno[[1]] |> filter(cluster=="Cluster 3"|
                                     cluster=="Cluster 4") |> 
  pull(Protein_id)

up <- c(t2v0_up[[1]],t4v0_up[[1]],
        t6v0_up[[1]],t8v0_up[[1]])
down <- c(t2v0_down[[1]],t4v0_down[[1]],
          t6v0_down[[1]],t8v0_down[[1]])

c12[!c12 %in% down]
#Not in down: "142158-jgi-mmetsp"

c34[!c34 %in% up]

up[!up %in% c34]
#Not in c34: "142158-jgi-mmetsp"

