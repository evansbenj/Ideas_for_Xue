library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(DESeq2) #differential expression analysis
library(tximport) #data import
library(readr) #read in data
library("RColorBrewer") #color for plot
library(gplots) #plotting

input_file_path <- file.path("D:", "school_grad", 
                             "Xue_master_projects/Project_borealis_sex_determination/Analysis",
                             "sample_clustering_analysis", "identify_by_borealisGonad_DEgene")

output_file_path <- input_file_path

input_file <- paste(input_file_path,"borealisGonadDEgene_tropTrna_blastnOut.tsv", sep="/")

input_data <- read_tsv(input_file, col_names = FALSE)


ortho_list <- input_data%>%
  group_by(X1)%>%
  arrange(X11, .by_group = TRUE)%>%
  top_n(1,X11)%>%
  top_n(1,X12)

rm(input_data)

ortho_trans <- ortho_list %>% 
  ungroup() %>% 
  select(X2) %>% 
  rename("Transcript_name"="X2")%>%
  unique()

expression_file = "../PCA/tropicalis_gonad.isoform.TMM.EXPR.matrix"

expression_data <- read.table(expression_file,header = TRUE, row.names ="Transcript_name")

ortho_expression <- expression_data %>%
  rownames_to_column("Transcript_name")%>%
  inner_join(ortho_trans, by=("Transcript_name"))%>%
  column_to_rownames("Transcript_name")

all_data <- DESeqDataSetFromMatrix(ortho_expression)
