#!/usr/bin/env Rscript
#----------------------------Usage-------------------------
#usage:  Rscript borealis_adult_de_R expression_matrix tissue_type output_folder
#
#-----------------------------------------------------------


#-------------------------------------getting arguments from command line-------------------------------------
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 2) {
  stop("At least two argument must be supplied (input file).n
       usage:  Rscript borealis_adult_de_R expression_matrix tissue_type output_folder", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  args[3] = pasete(args[1], "_de_result.tsv")
}

#--------------------------------------------------------------------------
#loading package
library(dplyr)
library(readr)
library(edgeR)

# define input and output files
expression_input = args[1]
de_output = args[2]
tissue_type = args[3]


#read in the data
express_data <- read.table(expression_input, header=T, row.names = 1, com='')

#sample names
#adult_liver_sample <- c("BJE3897_mom_liver","BJE4009_girl_liver","BJE4072_girl_liver","BJE4082_girl_liver","BJE3896_dad_liver","BJE3929_boy_liver","BJE4017_boy_liver","BJE4039_boy_liver")
liver_sample <- c("BJE3897_mom","BJE4009_girl","BJE4072_girl","BJE4082_girl","BJE3896_dad","BJE3929_boy","BJE4017_boy","BJE4039_boy")
gonad_sample <- c("BJE4009_girl_oviduct","BJE4072_girl_oviduct","BJE4082_girl_oviduct","BJE3896_dad_testis","BJE3929_boy_testis","BJE4017_boy_testis","BJE4039_boy_testis")
stage46 <- c("XBO8","XBO15","XBO24","XBO27", "XBO12","XBO16","XBO23","XBO26")
stage48 <- c("XBO19","XBO20","XBO29","XBO30","XBO33","XBO34","XBO35","XBO36","XBO17","XBO28","XBO31","XBO32")
stage48_49 <- c("XBO19","XBO20","XBO21","XBO29","XBO30","XBO33","XBO34","XBO35","XBO36","XBO17","XBO28","XBO31","XBO32")

#sex information
liver_sex <- factor(c(rep("female", 4), rep("male", 4)))
gonad_sex <- factor(c(rep("female", 3), rep("male", 4)))
stage46_sex <- factor(c(rep("female", 4), rep("male", 4)))
stage48_49_sex <- factor(c(rep("female", 8), rep("male", 5)))
stage48_sex <- factor(c(rep("female", 7), rep("male", 5)))

if (tissue_type=="liver"){
  sample = liver_sample
  conditions = liver_sex
  sample_size = 8
  female_size = 4
}

if(tissue_type=="gonad"){
  sample = gonad_sample
  conditions = gonad_sex
  sample_size = 7
  female_size = 3
}

if(tissue_type=="gonad_stage46"){
  sample = stage46
  conditions = stage46_sex
  sample_size = 8
  female_size = 4
}

if(tissue_type=="gonad_stage48"){
  sample = stage48_49
  conditions = stage48_49_sex
  sample_size = 13
  female_size = 8
}

#data filtering
express_data <- express_data[,match(sample, names(express_data))]
express_data <- express_data[rowSums(express_data)>= sample_size,] 
rnaseqMatrix <- express_data[!is.na(rowSums(express_data)),]

# do DE with edgeR
exp_study = DGEList(counts=rnaseqMatrix, group=conditions)
exp_study = calcNormFactors(exp_study)
exp_study = estimateCommonDisp(exp_study)
exp_study = estimateTagwiseDisp(exp_study)
et = exactTest(exp_study, pair=c("female", "male"))
tTags = topTags(et,n=NULL, sort.by = "none") # extract result
result_table = tTags$table
result_table = data.frame(sampleA="female", sampleB="male", result_table)
# exactTest: Note that the first group listed in the pair is the baseline for the 
# comparison—so if the pair is c("A","B") then the comparison is B - A, so 
# genes with positive log-fold change are up-regulated in group B compared 
# with group A (and vice versa for genes with negative log-fold change)

de_output <- paste0(output_folder,"de_edgeR_result_bor_",tissue_type,".tsv")

write.table(result_table, file =de_output, append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)



