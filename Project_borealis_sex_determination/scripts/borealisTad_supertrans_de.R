library(tidyverse)
library(graphics)
library(dplyr)
library(edgeR)


#list of input_files for borealis

input_file_folder <- file.path("Project_borealis_sex_determination/data/xb_transcriptome/supertrans/")

output_file_folder <- input_file_folder

input_file <- paste(input_file_folder,"borealisTad_supertrans_laevisGenome92_gmap_bedfile.bed", sep="/")

#read in the data
input_data <- read_tsv(input_file, col_names = FALSE)

#filter the mapping and select one mapping per transcripts
mapping <- (input_data 
            %>% select(chromosome = X1, 
                      start = X2, 
                       end = X3,
                       supertrans_id = X4,
                       mapping_quality = X5,
                       strandness = X6)
            %>% group_by(supertrans_id)
            %>% top_n(n=1, wt = mapping_quality)
            %>% filter (mapping_quality !=0)
            %>% mutate (count = count_)
)

#------------------------------de--------------------------------------------
#read in the data
input_file <- paste(input_file_folder,"borealis_tad_gonad_supertrans.counts.matrix", sep="/")
express_data <- read.table(input_file, header=T, row.names=1, com='')
express_data = express_data[rowSums(express_data)> 8,]


#stage information
stage46 <- c("XBO8","XBO15","XBO24","XBO27", "XBO12","XBO16","XBO23","XBO26")
stage48 <- c("XBO19","XBO20","XBO29","XBO30","XBO33","XBO34","XBO35","XBO36","XBO17","XBO28","XBO31","XBO32")
stage48_49 <- c("XBO19","XBO20","XBO21","XBO29","XBO30","XBO33","XBO34","XBO35","XBO36","XBO17","XBO28","XBO31","XBO32")

#sex information
stage46_sex <- factor(c(rep("female", 4), rep("male", 4)))
stage48_49_sex <- factor(c(rep("female", 8), rep("male", 5)))
stage48_sex <- factor(c(rep("female", 7), rep("male", 5)))

#select stage and conditions
rnaseqMatrix <- express_data[ ,which (names(express_data) %in% stage46)]
conditions = stage46_sex

# do DE with edgeR
exp_study = DGEList(counts=rnaseqMatrix, group=conditions)
exp_study = calcNormFactors(exp_study)
exp_study = estimateCommonDisp(exp_study)
exp_study = estimateTagwiseDisp(exp_study)
et = exactTest(exp_study, pair=c("female", "male"))
tTags = topTags(et,n=NULL)
result_table = tTags$table
result_table = data.frame(sampleA="female", sampleB="male", result_table)

#---------------------filter DE result-------------------------------------

summary (result_table)

result_table <- rownames_to_column(result_table, var = "supertrans_id")
result_filter <- (result_table
                  %>% filter(FDR <= 0.05)

) 

joined <- inner_join(result_filter, mapping, by = "supertrans_id")

#---------------------combine DE and Mapping-------------------------------------


