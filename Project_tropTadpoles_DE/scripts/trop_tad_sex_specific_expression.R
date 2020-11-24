library(tidyverse)
library(edgeR)

#read in the data

#transcriptome
input_file_folder <- file.path("project_sex_determination_gene/data/tropicalis_tad/transcriptome/raw_count_kallisto")
input_file <- paste(input_file_folder,"tropicalis_gonad.isoform.counts.matrix", sep="/")

express_data <- read.table(input_file, header=T, row.names = 1, com='')
#express_data = express_data[rowSums(express_data)> 10,]
express_data <- express_data[complete.cases(express_data), ]

#stage information
stage50 <- c("XT2","XT3","XT6","XT9", "XT10","XT11","XT16","XT17", "XT20","XT1","XT7","XT8","XT13","XT19" )


#sex information
stage50_sex <- factor(c(rep("female", 9), rep("male", 5)))

sex_specific <- (express_data
                 %>% rownames_to_column(var = "trans_id")
                 %>% mutate(
                   female_sum = XT2+XT3+XT6+XT9+XT10+XT11+XT16+XT17+XT20,
                   male_sum = XT1+XT7+XT8+XT13+XT19
                 )
                 %>%mutate( sex_specific = case_when(
                   female_sum == 0 & male_sum >0 ~ "male",
                   female_sum > 0 & male_sum ==0 ~ "female",
                   female_sum > 0 & male_sum >0 ~ "both",
                   female_sum == 0 & male_sum == 0 ~ "not_expressed"
                 )
                 )
                 %>% select(gene_id = trans_id, sex_specific, female_raw_sum = female_sum, male_raw_sum = male_sum)
)


output_folder <- file.path("project_sex_determination_gene/analysis/trop_tad/trop_tad_sex_specific_expression")
output_file <- paste(output_folder,"sex_specific_expression_list_trop_transcript.tsv", sep="/")
write.table(sex_specific, output_file, na = " ", sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)

