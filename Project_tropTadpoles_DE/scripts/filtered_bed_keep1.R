library(tidyverse)
library(graphics)
library(dplyr)
library(edgeR)

#-------------------------------------mapping information-------------------------------------

#list of input_files for borealis

#input_file_folder <- file.path("Project_borealis_sex_determination/data/xb_transcriptome/supertrans/")
#input_file <- paste(input_file_folder,"borealisTad_denovoT_laevisGenome92_gmap_bedfile.bed", sep="/")

input_file_folder <- file.path("Project_borealis_sex_determination/data/analysis/transcriptome/mapping_trans_laevisGenome92_gmap")
input_file <- paste(input_file_folder,"borealisTad_denovoT_laevisGenome92_gmap_bedfile.bed", sep="/")


output_folder <- input_file_folder

#read in the data
input_data <- read_tsv(input_file, col_names = FALSE)

#filter the mapping and select one mapping per transcripts
mapping <- (input_data 
            %>% select(chromosome = X1,  
                       start = X2, 
                       end = X3,
                       supertrans_id = X4,
                       mapping_quality = X5,
                       strandness = X6) # name the columns 
            %>% group_by(supertrans_id) 
            %>% top_n(n=1, wt = mapping_quality) # only select the one with the highest mapping quality
            %>% add_count(supertrans_id, name = "count")
            %>% mutate(chrCount = sum(!grepl("Scaffold", chromosome) ),
                       scaCount = sum(grepl("Scaffold", chromosome)))
)

SingleMatch <- (mapping 
                %>% filter(count == 1)
                %>% select(supertrans_id, chromosome, start, end)
)


mulMatch_sameChr <- (mapping
                     %>% group_by(supertrans_id)
                     %>% filter(count >1, count==chrCount, n_distinct(chromosome)==1)
                     %>% top_n(n=1, wt=start)
                     %>% select(supertrans_id, chromosome, start, end)
)


mulMatch_allSca <- (mapping
                    %>% filter(count >1, count == scaCount)
                    %>% group_by(supertrans_id)
                    %>% top_n(n=1, wt=chromosome)
                    %>% select(supertrans_id, chromosome, start, end)
)

mulMatch_OneChr <- (mapping
                    %>% filter(count >1, chrCount == 1, count > chrCount, 
                               !grepl("Scaffold", chromosome))
                    %>% select(supertrans_id, chromosome, start, end)
)

filtered_mapping <- bind_rows(SingleMatch,mulMatch_sameChr, mulMatch_allSca, mulMatch_OneChr)


output_file <- paste(output_folder,"borealisTad_denovoT_laevisV92_gmap_filtered.bed", sep="/")
write_tsv(filtered_mapping, output_file, na = "NA", append = FALSE, col_names = FALSE,
          quote_escape = "none")
