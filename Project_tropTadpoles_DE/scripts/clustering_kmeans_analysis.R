setwd("D:/school_grad school/Project_borealis_sex_determination/Analysis/sample_clustering_analysis/PCA/")
#setwd("~/Desktop/borealis_sex_determination/PCA")

library(datasets)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

#list of input_files
#input_file = "tropicalis_gonad.gene.counts.matrix"	
#input_file = "tropicalis_gonad.gene.TMM.EXPR.matrix"
#input_file = "tropicalis_gonad.isoform.counts.matrix"
input_file = "tropicalis_gonad.isoform.TMM.EXPR.matrix"

#read in the data
input_data <- read.table(input_file,header = TRUE, row.names =1)

express_data <- input_data[rowSums(input_data)>0,]
express_data <- express_data[,]
express_data <- t(express_data)
express_data <- scale(express_data)

set.seed(20)
expCluster <- kmeans(express_data, 2, nstart = 20)
expCluster$cluster
expCluster$iter
expCluster$ifault


fviz_cluster(expCluster, data = express_data)

