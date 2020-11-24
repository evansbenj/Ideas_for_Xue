library(tidyverse)
library(graphics)
library(dplyr)
library(edgeR)
library(IRanges)
library(GenomicRanges)

#-------------------------------------collapsing expression levels-------------------------------------
#read in mapping data
input_file_folder <- file.path("Project_borealis_sexual_antagonism/data/transcriptome/transcriptome_dec2018/mapping_borealisTrans_borealisGenome_gmap")
input_file <- paste(input_file_folder,"borealis_denovoT_borealis_genome_gmap_collapsed.tsv", sep="/")
collaspse_transcript_list <- (read_tsv(input_file, col_names = TRUE) 
                              %>% mutate_at(vars(matches("bin_id")),list(factor))
)

#read in expression data
input_file_folder <- file.path("Project_borealis_sexual_antagonism/data/transcriptome/transcriptome_dec2018/raw_count_salmon")
input_file <- paste(input_file_folder,"borealis_adult_liver.counts.matrix", sep="/")
expression <- read_tsv(input_file, col_names = TRUE) 

#specific output file
output_folder <- input_file_folder

#collapse expression value for each transcript bin
collapse_expression <- (left_join(collaspse_transcript_list, expression, by =c("trans_id"="X1"))
                        %>% select(-trans_id)
                        %>% group_by(bin_id)
                        %>% summarise_each(list(sum))
)

output_file <- paste(output_folder,"borealis_adult_liver_collapsedCounts.matrix", sep="/")
write_tsv(collapse_trans_list, output_file, na = "NA", append = FALSE, col_names = TRUE,
          quote_escape = "none")


                     