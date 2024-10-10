# Download and prepare data from TCGA
# Joel Ruiz Hernandez
# Contact: jjrhernandez@gmail.com
# Description: This script downloads RNA-Seq data from the TCGA-CESC project and create metadata.
# Output: A count matrix (TCGA_CESC_Data_matrix) saved as a TSV file.

## Set working directory
setwd("~/1_Get_Data_TCGA/")

#load("1_Image_data_RNA.RData")

## Packages required
library(TCGAbiolinks)
library(tidyverse)
library(SummarizedExperiment)
library(vroom)

## Get query (met tumor excluded)
sample_types <- c("Primary Tumor", "Solid Tissue Normal")

query_TCGA_CESC <- GDCquery(
  project = "TCGA-CESC",
  data.category = 'Transcriptome Profiling',
  experimental.strategy = 'RNA-Seq',
  workflow.type = 'STAR - Counts',
  access = 'open',
  sample.type = sample_types)

Output_query_TCGA <- getResults(query_TCGA_CESC) %>% 
  dplyr::select(cases, cases.submitter_id, sample_type)
#head(Output_query_TCGA)
#   cases                              cases.submitter_id  sample_type
#1 TCGA-C5-A1M5-01A-11R-A13Y-07       TCGA-C5-A1M5 Primary Tumor
#2 TCGA-EK-A2R9-01A-11R-A18M-07       TCGA-EK-A2R9 Primary Tumor


## Download data 
GDCdownload(query_TCGA_CESC) #307 files in GDC directory

## Prepare data
TCGA_prepare <- GDCprepare(query_TCGA_CESC, summarizedExperiment = T)
    #Available assays in SummarizedExperiment : 
    #=> unstranded
    #=> stranded_first
    #=> stranded_second
    #=> tpm_unstrand
    #=> fpkm_unstrand
    #=> fpkm_uq_unstrand
unstranded <- assay(TCGA_prepare, 'unstranded')
unstranded <- as_tibble(unstranded, rownames = "gene")
#dim(unstranded) #[1] 60660   308
#head(unstranded[1:10,1:3])
#gene               `TCGA-C5-A1M5-01A-11R-A13Y-07` `TCGA-EK-A2R9-01A-11R-A18M-07`
#<chr>                                       <int>                          <int>
#1 ENSG00000000003.15                           6914                          10140
#2 ENSG00000000005.6                               5                              1
#3 ENSG00000000419.13                           3793                           4078



##Create metadata --- ---

#Add HPV type 
HPV_TCGA <- vroom(file = "~/0_HPV_Distribution/1_2_HPV_Clades.tsv")  #[1] 275   3
#head(HPV_TCGA)
#`TCGA Case ID` tipos_HPV Clado_filogenetico
#<chr>          <chr>     <chr>             
#1 TCGA-HM-A4S6   HPV16     A9                
#2 TCGA-ZJ-AAXD   HPV16     A9                
#3 TCGA-EA-A43B   HPV16     A9 

#CLinical info
Clinical_info <- GDCquery_clinic("TCGA-CESC","clinical") %>% 
  dplyr::select("submitter_id", "figo_stage", "race","primary_diagnosis")
# Bind
Metadata <- Output_query_TCGA %>% 
  left_join(Clinical_info, by = c("cases.submitter_id" = "submitter_id")) %>% 
  right_join(HPV_TCGA, by = c("cases.submitter_id" = "TCGA Case ID")) %>%
  mutate_at(vars(figo_stage:Clado_filogenetico, -race), 
            ~ ifelse(sample_type == "Solid Tissue Normal", "Solid Tissue Normal", .)) %>% 
  dplyr::rename(specimenID =  cases) %>% dplyr::rename(HPV_type = tipos_HPV) %>% 
  dplyr::rename(HPV_clade = Clado_filogenetico)
#[1] 278   8
#colnames(Metadata)
#[1] "specimenID"         "cases.submitter_id" "sample_type"        "figo_stage"         "race"              
#[6] "primary_diagnosis"  "HPV_type"           "HPV_clade" 

#Filter unstranded with the cases in Metadata
Cases_metadata <- Metadata %>% pull(specimenID)
unstranded_counts <- unstranded %>% dplyr::select(gene, Cases_metadata)
    #dim(unstranded)#[1] 60660   308
    #dim(unstranded_counts) #[1] 60660   279   #Note: The 29 left were HPV negative samples
    #verify: all(Cases_metadata %in% colnames(unstranded))  #TRUE

# Num samples by HPV_clade
CladeHPV_specimenID <- Metadata %>%
  group_by(HPV_clade) %>%
  summarise(num_muestras = n()) %>%
  pivot_wider(names_from = HPV_clade, values_from = num_muestras)
#A7    A9 negative  otro
#<int> <int>    <int> <int>
#  1    68   202        3     5

##Save metadata
Metadata <- as_tibble(Metadata)
vroom_write(Metadata, "1_2_Metadata.tsv", delim = "\t")

## Save count matrix
unstranded_counts <- as_tibble(unstranded_counts)
vroom_write(unstranded_counts, "1_1_unstranded_counts.tsv", delim = "\t")

#Save work space
save.image(file = "1_Image_data_RNA.RData")

