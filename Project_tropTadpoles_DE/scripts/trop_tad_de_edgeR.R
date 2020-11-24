library(tidyverse)
library(edgeR)

#read in the data

#transcriptome
input_file_folder <- file.path("project_sex_determination_gene/data/tropicalis_tad/raw_count_kallisto")
input_file <- paste(input_file_folder,"tropicalis_gonad.isoform.counts.matrix", sep="/")

express_data <- read.table(input_file, header=T, row.names = 1, com='')
express_data <- express_data[complete.cases(express_data), ]

#stage information
stage50 <- c("XT2","XT3","XT6","XT9", "XT10","XT11","XT16","XT17", "XT20","XT1","XT7","XT8","XT13","XT19" )
group1 <- c("XT3", "XT9", "XT11","XT20", "XT7", "XT8")
group2 <- c("XT2", "XT6", "XT10", "XT16", "XT17", "XT1", "XT13", "XT19")
group3 <- c("XT3", "XT9", "XT11","XT20","XT1", "XT13", "XT19")
group4 <- c("XT2", "XT6","XT10", "XT16", "XT17","XT7", "XT8")
#sex information
stage50_sex <- factor(c(rep("female", 9), rep("male", 5)))
group1_sex <- factor(c(rep("female", 4), rep("male", 2)))
group2_sex <- factor(c(rep("female", 5), rep("male", 3)))
group3_sex <- factor(c(rep("female", 4), rep("male", 3)))
group4_sex <- factor(c(rep("female", 5), rep("male", 2)))

length(stage)

#select stage and conditions
stage = group4
conditions = group4_sex
rnaseqMatrix <- express_data[,match(stage, names(express_data))]
rnaseqMatrix = rnaseqMatrix[rowSums(rnaseqMatrix)> length(stage),]


# do DE with edgeR
exp_study = DGEList(counts=rnaseqMatrix, group=conditions)
exp_study = calcNormFactors(exp_study)
exp_study = estimateCommonDisp(exp_study)
exp_study = estimateTagwiseDisp(exp_study)
et = exactTest(exp_study, pair=c("female", "male"))
tTags = topTags(et,n=NULL, sort.by = "none")
result_table = tTags$table
result_table = data.frame(sampleA="female", sampleB="male", result_table)

#exp_study$samples$group


output_folder <- file.path("project_sex_determination_gene/data/tropicalis_tad/trop_tad_de_result")
output_file <- paste(output_folder,"de_result_trop_uncollapsed_group4_edgeR.tsv", sep="/")
write.table(result_table, output_file, na = " ", sep = "\t", quote = FALSE, row.names = TRUE,col.names = TRUE)



