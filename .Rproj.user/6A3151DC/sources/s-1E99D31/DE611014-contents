#---------------------------------------perform DE--------------------------
library("limma")
library("edgeR")
library(tidyverse)

#list of input_files for borealis
input_file = "Project_borealis_sex_determination/data/xb_transcriptome/counts_kallisto/borealis_tad_gonad_supertrans.counts.matrix"
#input_file = "Project_borealis_sex_determination/data/xb_transcriptome/counts_kallisto/borealis_tad_gonad_trans.counts.matrix"


#read in the data
input_data <- read.table(input_file,header = TRUE)

#do some filtering
#filtering method 1
express_data = input_data[rowSums(input_data)> 1,]
#filtering method 2
#express_data = input_data[rowSums(cpm(input_data) > 1) >= 1,]
#filtering method 3
# cpm_log <- cpm(rnaseqMatrix, log = TRUE)
# median_log2_cpm <- apply(cpm_log, 1, median)
# expr_cutoff <- -3
# rnaseqMatrix <- rnaseqMatrix[median_log2_cpm > expr_cutoff, ]

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

#---------------------examine result-------------------------------------

summary (result_table)

result_table <- rownames_to_column(result_table)
result_filter <- (result_table
                  %>% filter(FDR <= 0.05)
)


