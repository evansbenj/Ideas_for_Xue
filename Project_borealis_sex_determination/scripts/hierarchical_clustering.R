library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(DESeq2) #differential expression analysis
library(tximport) #data import
library(readr) #read in data
library("RColorBrewer") #color for plot
library(gplots) #plotting

input_file_path <- file.path("D:", "school_grad school", 
                             "Xue_master_projects/Project_borealis_sex_determination/Analysis",
                             "sample_clustering_analysis", "identify_by_borealisGonad_DEgene")

output_file_path <- input_file_path

