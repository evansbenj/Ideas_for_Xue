#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = pasete(args[1], ".keep1match.bed")
}



#-------------------------------------mapping information-------------------------------------
library(dplyr)
library(readr)

# define input and output files
mapping_input = args[1]
mapping_output = args[2]

#read in the data
input_data <- read_tsv(mapping_input, col_names = FALSE)

#filter the mapping and select one mapping per transcripts
mapping <- (input_data 
            %>% select(chromosome = X1,  
                       start = X2, 
                       end = X3,
                       supertrans_id = X4,
                       mapping_quality = X5,
                       strandness = X6) # name the columns 
            %>% group_by(supertrans_id) 
            %>% top_n(n=1, wt = mapping_quality) # only select the one with the highest mapping quality; would keep ties
            %>% add_count(supertrans_id, name = "count")
            %>% mutate(chrCount = sum(!grepl("ca", chromosome) ), #scaffold can start with S or s; omit the first letter to avoid confusion
                       scaCount = sum(grepl("ca", chromosome)))
)

#str(mapping)

SingleMatch <- (mapping 
                %>% filter(count == 1)
                %>% select(trans_id = supertrans_id, chromosome, start, end)
)

#str(SingleMatch)

mulMatch_sameChr <- (mapping
                     %>% group_by(supertrans_id)
                     %>% filter(count >1 & count==chrCount & n_distinct(chromosome)==1)
                     %>% top_n(n=1, wt=start)
                     %>% select(trans_id = supertrans_id, chromosome, start, end)
)

#str(mulMatch_sameChr)

mulMatch_allSca <- (mapping
                    %>% filter(count >1 & count == scaCount)
                    %>% group_by(supertrans_id)
                    %>% top_n(n=1, wt=chromosome)
                    %>% top_n(n=1, wt=start)
                    %>% select(trans_id = supertrans_id, chromosome, start, end)
)

#str(mulMatch_allSca)

mulMatch_OneChr <- (mapping
                    %>% filter(count >1 & chrCount == 1 & count > chrCount & 
                               !grepl("ca", chromosome))
                    %>% select(trans_id = supertrans_id, chromosome, start, end)
)

#str(mulMatch_allSca)

filtered_mapping <- bind_rows(SingleMatch,mulMatch_sameChr, mulMatch_allSca, mulMatch_OneChr)

write_tsv(filtered_mapping, mapping_output, na = "NA", append = FALSE, col_names = FALSE,
          quote_escape = "none")
