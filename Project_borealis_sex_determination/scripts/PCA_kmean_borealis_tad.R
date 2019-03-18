library(tidyverse)
library(graphics)

#another method Caroline recommended 
#dudi.pca
#library(mass)
#library(ad4)

#list of input_files for trop
#input_file = "tropicalis_gonad.gene.counts.matrix"	
#input_file = "tropicalis_gonad.gene.TMM.EXPR.matrix"
#input_file = "tropicalis_gonad.isoform.counts.matrix"
#input_file = "tropicalis_gonad.isoform.TMM.EXPR.matrix"

#list of input_files for borealis
input_file = "Project_borealis_sex_determination/data/borealis_tad_gonad_supertrans.TMM.EXPR.matrix"
input_file = "Project_borealis_sex_determination/data/borealis_tad_gonad_trans.TMM.EXPR.matrix"


#read in the data
input_data <- read.table(input_file,header = TRUE, row.names ="Transcript_name")

#do some filtering
# I filtered out the genes with zeros expression in all samples;
express_data <- input_data[rowSums(input_data)>1,]
#stage information
stage46 <- c("XBO12","XBO15","XBO16","XBO23","XBO24","XBO26","XBO27","XBO08")
stage48 <- c("XBO17","XBO19","XBO20","XBO21","XBO28","XBO29","XBO30","XBO31","XBO32","XBO33","XBO34","XBO35","XBO36")

express_data <- express_data[,which (names(express_data) %in% stage46)]

#prcomp would expect the column to be variable (different genes) and row to be samples. Below is an example.
#Sample Gene1 Gene2 Gene3 Gene4
#Sample1  5 5 10`10`
#Sample2  2 3 15  29
#Hence I transcponse the data into the format that prcomp expects
express_data <- t(express_data)
express_data <- express_data[,colSums(express_data)>0]

#Do PCA analysis using prcomp 
#By default, center = TRUE and it centers the variable to have mean equals to zero. 
#With parameter scale. = TRUE, we normalize the variables to have standard deviation equals to 1.
#prcomp return the following items:
###sdev: the standard deviations of the principal components
###rotation: the matrix of variable loadings (i.e., a matrix whose columns contain the eigenvectors). 
###x: if retx is true the value of the rotated data (the centred (and scaled if requested) data multiplied by the rotation matrix) is returned
###center, scale: the centering and scaling used
project_pca <- prcomp(express_data, retx=TRUE, center = TRUE, scale. = TRUE)

#calculate the project_pca_proportionvarianes
project_pca_proportionvariances <- ((project_pca$sdev^2) / (sum(project_pca$sdev^2)))*100

#look at the pca result
summary(project_pca)

################################Visualization############################################
## now make a fancy looking plot that shows the PCs and the variation:
library(ggplot2)

pca_plot <- data.frame(Sample=rownames(project_pca$x),
                       X=project_pca$x[,1],
                       Y=project_pca$x[,2])
pca_plot

pca_plot<-ggplot(data=pca_plot, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", project_pca_proportionvariances[1], "%", sep="")) +
  ylab(paste("PC2 - ", project_pca_proportionvariances[2], "%", sep="")) +
  theme_classic() +
  ggtitle("Gonad RNAseq PCA") +
  geom_hline(mapping = NULL, data = NULL, yintercept = 0,
             na.rm = FALSE, show.legend = NA) +
  geom_vline(mapping = NULL, data = NULL, xintercept = 0,
             na.rm = FALSE, show.legend = NA)

pca_plot

express_data <- scale(express_data)

set.seed(20)
expCluster <- kmeans(express_data, 2, nstart = 20)
expCluster$cluster
expCluster$iter
expCluster$ifault


kmean_plot <- fviz_cluster(expCluster, data = express_data)
kmean_plot
