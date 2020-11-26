#!/usr/bin/env Rscript
#---------------------Usage----------------------------------------------------------
#
# usage:  Rscript bor_de_deseq2.R path_to_expression_matrix tissue_type output_folder
#
#--------------------getting arguments from command line-------------------------------------

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 2) {
  stop("At least two argument must be supplied (input file).n
       usage:  Rscript bor_de_deseq2.R path_to_expression_matrix tissue_type output_folder", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  args[3] = pasete(args[1], "_de_result.tsv")
}

#--------------------load package-----------------------------------------------------
library(limma)
library(edgeR)
library(DESeq2)
library(apeglm)


# define input and output files: from user
expression_input = args[1]
tissue_type = args[2]
output_folder = args[3]

express_data <- read.table(expression_input, header=T, row.names = 1, com='')

#stage information
liver_sample <- c("BJE3897_mom_liver","BJE4009_girl_liver","BJE4072_girl_liver","BJE4082_girl_liver","BJE3896_dad_liver","BJE3929_boy_liver","BJE4017_boy_liver","BJE4039_boy_liver")
gonad_sample <- c("BJE4009_girl_oviduct","BJE4072_girl_oviduct","BJE4082_girl_oviduct","BJE3896_dad_testis","BJE3929_boy_testis","BJE4017_boy_testis","BJE4039_boy_testis")
stage46 <- c("XBO8","XBO15","XBO24","XBO27", "XBO12","XBO16","XBO23","XBO26")
stage48 <- c("XBO19","XBO20","XBO29","XBO30","XBO33","XBO34","XBO35","XBO36","XBO17","XBO28","XBO31","XBO32")
stage48_49 <- c("XBO19","XBO20","XBO21","XBO29","XBO30","XBO33","XBO34","XBO35","XBO36","XBO17","XBO28","XBO31","XBO32")

#sex information
liver_sex <- data.frame(conditions=factor(c(rep("female", 4), rep("male", 4))))
gonad_sex <- data.frame(conditions=factor(c(rep("female", 3), rep("male", 4))))
stage46_sex <- data.frame(conditions=factor(c(rep("female", 4), rep("male", 4))))
stage48_49_sex <- data.frame(conditions=factor(c(rep("female", 8), rep("male", 5))))
stage48_sex <- data.frame(conditions=factor(c(rep("female", 7), rep("male", 5))))

if (tissue_type=="liver"){
  sample = liver_sample
  conditions = liver_sex
  sample_size = 8
}

if(tissue_type=="gonad"){
  sample = gonad_sample
  conditions = gonad_sex
  sample_size = 7
}

if(tissue_type=="gonad_stage46"){
  sample = stage46
  conditions = stage46_sex
  sample_size = 8
}

if(tissue_type=="gonad_stage48"){
  sample = stage48_49
  conditions = stage48_49_sex
  sample_size = 13
}

#filter and prepare for deseq object 
express_data <- express_data[complete.cases(express_data), ]
rnaseqMatrix <- express_data[ ,match(sample,names(express_data))]
rnaseqMatrix <- rnaseqMatrix[rowSums(rnaseqMatrix)>sample_size,]
rnaseqMatrix = round(rnaseqMatrix) #DESeq2 expect integers
rownames(conditions) = colnames(rnaseqMatrix)


#store it in the DESeq DE object
#design: y ~ x1+x2+x3+sex, x1 are the batch effect like if they are in different lane
#biological variable
#2 sexes: 4 individual for sex 
#include/exclude parental samples
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = rnaseqMatrix,
  colData = conditions, #a table with metadata on the count tables columns
  design = ~ conditions) #design formula is used to estimate the dispersions and to estimate the log2 fold changes of the model

#the function DESeq, what it does:
#1)estimation of size factors (which control for differences 
#in the library size of the sequencing experiments),
#2) the estimation of dispersion for each gene, 
#3)and fitting a generalized linear model.

dds = DESeq(ddsFullCountTable)
#plotDispEsts((dds))

#contrast:t specifies what comparison to extract from the object to build a results table
#simplest case: a character vector with exactly three elements: 
#1. the name of a factor in the design formula, -----> condition
#2. the name of the numerator level for the fold change, and ---> female
#3. the name of the denominator level for the fold change --->male
#we don't need contrast because we only have one comparison: male vs female
#contrast=c("conditions","female","male")

#Extract results from a DESeq analysis
# res <- results(dds, contrast = contrast)
# res #condition: male vs female; it means the log2FC is log2(male/female)
# summary(res)
# res <- res[order(res$padj),]


resLFC <- lfcShrink(dds, coef = 2, type = "apeglm")
summary(resLFC)


de_output <- paste0(output_folder,"de_deseq2_result_bor_uncollapsed_",tissue_type,".tsv")
write.table(resLFC, file = de_output, append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)

##### plot
# ma <- plotMA(resLFC, ylim =c(-5,5), cex = 0.5)
# ma
# 
# head(resLFC, n = 20)
# resLFC$padj[is.na(resLFC$padj)]  <- 1 #fill the NA with 1 so that it will later filter out
# res_de_only= resLFC[resLFC$padj < 0.05,]
# #res = res[res$log2FoldChange >= 1|res$log2FoldChange <= -1,]
# plotMA(resLFC, ylim =c(-5,5), cex = 0.5)