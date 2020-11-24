library(tidyverse)
library(DESeq2)
library(apeglm)

#read in the data

#count input
input_file_folder <- file.path("project_sex_determination_gene/data/tropicalis_tad/transcriptome/raw_count_kallisto")
input_file <- paste(input_file_folder,"tropicalis_gonad.isoform.counts.matrix", sep="/")

express_data <- read.table(input_file, header=T, row.names = 1, com='')

##sample information
#stage information
stage50 <- c("XT2","XT3","XT6","XT9", "XT10","XT11","XT16","XT17", "XT20","XT1","XT7","XT8","XT13","XT19" )
#sex information
stage50_sex <- data.frame(conditions=factor(c(rep("female", 9), rep("male", 5))))
#sample size
sample_size = 14

#filter and prepare for deseq object 
express_data <- express_data[complete.cases(express_data), ]
express_data <- express_data[rowSums(express_data)>5,]
express_data <- express_data[!is.na(rowSums(express_data)),]
rnaseqMatrix <- express_data[,match(stage50, names(express_data))]
rnaseqMatrix = round(rnaseqMatrix) #DESeq2 expect integers
conditions = stage50_sex
rownames(conditions) = colnames(rnaseqMatrix)

#store it in the DESeq DE object
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = rnaseqMatrix,
  colData = conditions, 
  design = ~ conditions) 

#the function DESeq, what it does:
#1)estimation of size factors (which control for differences 
#in the library size of the sequencing experiments),
#2) the estimation of dispersion for each gene, 
#3)and fitting a generalized linear model.
dds = DESeq(ddsFullCountTable)
plotDispEsts((dds))

#shrinkage
resLFC <- lfcShrink(dds, coef = 2, type = "apeglm")
summary(resLFC)

#output
output_folder <- file.path("project_sex_determination_gene/analysis/trop_tad/trop_tad_de_result/")
de_output <- file.path(output_folder,"de_result_trop_stage50_deseq2.tsv")
write.table(resLFC, file = de_output, append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)

