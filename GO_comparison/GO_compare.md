GO\_Compare
================
Alex Gil
8/13/2021

# Comparing GO terms from an old and a new analysis

file location: /gpfs/projects/RestGroup/mariana/GO/

## The old analysis was generated with JGI

  - GO ID:
      - Aurli1\_GeneCatalog\_proteins\_20120618\_GO.tab

## The new analysis was generated with [sanspanz](http://ekhidna2.biocenter.helsinki.fi/sanspanz/)

  - GO ID:
      - GO\_predi\_det\_Aurli1\_ekhidna2.txt

## Load files

``` r
library("tidyverse")
```

    ## Warning: package 'tidyverse' was built under R version 4.0.5

    ## -- Attaching packages --------------------------------------- tidyverse 1.3.1 --

    ## v ggplot2 3.3.5     v purrr   0.3.4
    ## v tibble  3.1.3     v dplyr   1.0.7
    ## v tidyr   1.1.3     v stringr 1.4.0
    ## v readr   2.0.1     v forcats 0.5.1

    ## Warning: package 'ggplot2' was built under R version 4.0.5

    ## Warning: package 'tibble' was built under R version 4.0.5

    ## Warning: package 'readr' was built under R version 4.0.5

    ## Warning: package 'dplyr' was built under R version 4.0.5

    ## Warning: package 'forcats' was built under R version 4.0.4

    ## -- Conflicts ------------------------------------------ tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library("here")
```

    ## Warning: package 'here' was built under R version 4.0.5

    ## here() starts at C:/Users/Alex/GitHub/laby

``` r
library("knitr")
```

    ## Warning: package 'knitr' was built under R version 4.0.5

``` r
jgi <- read_delim(here("GO_comparison","Aurli1_GeneCatalog_proteins_20120618_GO.tab"),delim = "\t")
```

    ## Rows: 31047 Columns: 5

    ## -- Column specification --------------------------------------------------------
    ## Delimiter: "\t"
    ## chr (3): goName, gotermType, goAcc
    ## dbl (2): #proteinId, gotermId

    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
jgi_GO <- jgi %>% 
  select(goAcc) %>%
  unique() %>% 
  pull()

sans <- read_delim(here("GO_comparison","GO_predi_det_Aurli1_ekhidna2.txt"),delim = "\t")%>% rename(goAcc=goid) %>%
  mutate(goAcc=paste0("GO:",goAcc))
```

    ## Rows: 42580 Columns: 7

    ## -- Column specification --------------------------------------------------------
    ## Delimiter: "\t"
    ## chr (4): qpid, ontology, goid, desc
    ## dbl (3): ARGOT_score, ARGOT_PPV, ARGOT_rank

    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
sans_GO <- sans  %>%
  select(goAcc) %>% 
  unique() %>% 
  pull()
```

# How many unique GO terms were found in each analysis?

``` r
print("JGI")
```

    ## [1] "JGI"

``` r
jgi_GO  %>% length() #2156 unique Go_ids
```

    ## [1] 2156

``` r
print("Sanzpanz")
```

    ## [1] "Sanzpanz"

``` r
sans_GO  %>% length() #7234  unique Go_ids
```

    ## [1] 7234

# Which GO terms are in each set?

``` r
tot <- c(sans_GO,jgi_GO) #8060 unique GO terms in total
first <- setdiff(jgi_GO,sans_GO) #826 unique GO terms in JGI
second <- setdiff(sans_GO,jgi_GO)  #5904 unique GO terms in Sans
both <- intersect(sans_GO,jgi_GO) #1330 unique GO terms in both
```

## Genes that were present in analysis 1, but not in 2.

``` r
set1 <- jgi %>% filter(goAcc %in% first) %>% 
  select(goAcc,goName) %>%
  distinct() %>% 
  arrange(goAcc) 

set1
```

    ## # A tibble: 826 x 2
    ##    goAcc      goName                                               
    ##    <chr>      <chr>                                                
    ##  1 GO:0000059 protein import into nucleus, docking                 
    ##  2 GO:0000104 succinate dehydrogenase activity                     
    ##  3 GO:0000119 mediator complex                                     
    ##  4 GO:0000140 acylglycerone-phosphate reductase activity           
    ##  5 GO:0000156 two-component response regulator activity            
    ##  6 GO:0000214 tRNA-intron endonuclease complex                     
    ##  7 GO:0000221 vacuolar proton-transporting V-type ATPase, V1 domain
    ##  8 GO:0000247 C-8 sterol isomerase activity                        
    ##  9 GO:0000253 3-keto sterol reductase activity                     
    ## 10 GO:0000254 C-4 methylsterol oxidase activity                    
    ## # ... with 816 more rows

## Genes that were present in analysis 2, but not in 1.

``` r
set2 <- sans %>% filter(goAcc %in% second) %>% 
  select(goAcc,desc) %>%
  distinct() %>% 
  arrange(goAcc)

set2
```

    ## # A tibble: 5,904 x 2
    ##    goAcc      desc                               
    ##    <chr>      <chr>                              
    ##  1 GO:0000001 mitochondrion inheritance          
    ##  2 GO:0000002 mitochondrial genome maintenance   
    ##  3 GO:0000011 vacuole inheritance                
    ##  4 GO:0000018 regulation of DNA recombination    
    ##  5 GO:0000019 regulation of mitotic recombination
    ##  6 GO:0000022 mitotic spindle elongation         
    ##  7 GO:0000023 maltose metabolic process          
    ##  8 GO:0000027 ribosomal large subunit assembly   
    ##  9 GO:0000028 ribosomal small subunit assembly   
    ## 10 GO:0000035 acyl binding                       
    ## # ... with 5,894 more rows

## Genes that were present in both analyses.

``` r
set_both <- sans %>% filter(goAcc %in% both) %>% 
  select(goAcc,desc) %>%
  distinct() %>% 
  arrange(goAcc)

set_both
```

    ## # A tibble: 1,330 x 2
    ##    goAcc      desc                                  
    ##    <chr>      <chr>                                 
    ##  1 GO:0000009 alpha-1,6-mannosyltransferase activity
    ##  2 GO:0000015 phosphopyruvate hydratase complex     
    ##  3 GO:0000026 alpha-1,2-mannosyltransferase activity
    ##  4 GO:0000030 mannosyltransferase activity          
    ##  5 GO:0000033 alpha-1,3-mannosyltransferase activity
    ##  6 GO:0000036 acyl carrier activity                 
    ##  7 GO:0000049 tRNA binding                          
    ##  8 GO:0000062 fatty-acyl-CoA binding                
    ##  9 GO:0000103 sulfate assimilation                  
    ## 10 GO:0000105 histidine biosynthetic process        
    ## # ... with 1,320 more rows

## Save results.

``` r
write_csv(set1,here("GO_comparison","set1.csv"))
write_csv(set2,here("GO_comparison","set2.csv"))
write_csv(set_both,here("GO_comparison","set_both.csv"))
```

## Annotations

  - Is there a difference in the \# of genes that got annotations?

Yes, JGI got 6935 unique protein IDs, SansPanz got 5755 unique protein
IDs.

``` r
jgi_anno <- jgi %>% 
  rename(Protein_ID = `#proteinId`) %>% 
  select(Protein_ID) %>% 
  distinct() %>% 
  arrange(Protein_ID) %>% 
  pull() 
  

jgi_anno %>% length()  # 6935 protein IDs
```

    ## [1] 6935

``` r
sans_anno <- sans  %>%
  select(qpid) %>% 
  mutate(anno = str_split(qpid, "\\|",4)) %>% 
  mutate(Protein_ID = map(anno,nth,3)) %>% 
  mutate(Protein_ID=as.numeric(Protein_ID)) %>% 
  select(Protein_ID) %>% 
  distinct() %>% 
  arrange(Protein_ID) %>% 
  pull() 

sans_anno %>% length() # 5755 protein IDs
```

    ## [1] 5755

# Which annotations are in each set?

``` r
tot <- c(sans_anno,jgi_anno) #8567 unique protein  IDs in total
first <- setdiff(jgi_anno,sans_anno) #2812 unique GO terms in JGI
second <- setdiff(sans_anno,jgi_anno)  #1632 unique GO terms in Sans
both <- intersect(sans_anno,jgi_anno) #4123 unique GO terms in both
```

  - How many genes were annotated in the old and new versions?

4123

``` r
set_both <- sans  %>%
  select(qpid) %>% 
  mutate(anno = str_split(qpid, "\\|",4)) %>% 
  mutate(Protein_ID = map(anno,nth,3)) %>% 
  mutate(Protein_ID=as.numeric(Protein_ID)) %>% 
  select(Protein_ID) %>% 
  filter(Protein_ID %in% both) %>% 
  select(Protein_ID) %>%
  distinct() %>% 
  arrange(Protein_ID) %>% pull()
```

  - How many not annotated in the old are now annotated in the new?

1632

``` r
set2 <- sans  %>%
  select(qpid) %>% 
  mutate(anno = str_split(qpid, "\\|",4)) %>% 
  mutate(Protein_ID = map(anno,nth,3)) %>% 
  mutate(Protein_ID=as.numeric(Protein_ID)) %>% 
  select(Protein_ID) %>%
  filter(Protein_ID %in% second) %>% 
  select(Protein_ID) %>%
  distinct() %>% 
  arrange(Protein_ID) %>% pull()
```

  - How many not annotated in the new were annotated in the old?

2812

``` r
set1 <- jgi %>%  
  rename(Protein_ID = `#proteinId`) %>% 
  distinct() %>%
  filter(Protein_ID %in% first) %>% 
  select(goAcc) %>% 
  arrange(goAcc) %>% 
  pull()
```

``` r
sessionInfo()
```

    ## R version 4.0.3 (2020-10-10)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19043)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.1252 
    ## [2] LC_CTYPE=English_United States.1252   
    ## [3] LC_MONETARY=English_United States.1252
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.1252    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] knitr_1.33      here_1.0.1      forcats_0.5.1   stringr_1.4.0  
    ##  [5] dplyr_1.0.7     purrr_0.3.4     readr_2.0.1     tidyr_1.1.3    
    ##  [9] tibble_3.1.3    ggplot2_3.3.5   tidyverse_1.3.1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.1.1  xfun_0.25         haven_2.4.3       colorspace_2.0-2 
    ##  [5] vctrs_0.3.8       generics_0.1.0    htmltools_0.5.1.1 yaml_2.2.1       
    ##  [9] utf8_1.2.2        rlang_0.4.11      pillar_1.6.2      glue_1.4.2       
    ## [13] withr_2.4.2       DBI_1.1.1         bit64_4.0.5       dbplyr_2.1.1     
    ## [17] modelr_0.1.8      readxl_1.3.1      lifecycle_1.0.0   munsell_0.5.0    
    ## [21] gtable_0.3.0      cellranger_1.1.0  rvest_1.0.1       evaluate_0.14    
    ## [25] tzdb_0.1.2        parallel_4.0.3    fansi_0.5.0       broom_0.7.9      
    ## [29] Rcpp_1.0.7        scales_1.1.1      backports_1.2.1   vroom_1.5.4      
    ## [33] jsonlite_1.7.2    bit_4.0.4         fs_1.5.0          hms_1.1.0        
    ## [37] digest_0.6.27     stringi_1.7.3     rprojroot_2.0.2   grid_4.0.3       
    ## [41] cli_3.0.1         tools_4.0.3       magrittr_2.0.1    crayon_1.4.1     
    ## [45] pkgconfig_2.0.3   ellipsis_0.3.2    xml2_1.3.2        reprex_2.0.1     
    ## [49] lubridate_1.7.10  assertthat_0.2.1  rmarkdown_2.10    httr_1.4.2       
    ## [53] rstudioapi_0.13   R6_2.5.1          compiler_4.0.3
