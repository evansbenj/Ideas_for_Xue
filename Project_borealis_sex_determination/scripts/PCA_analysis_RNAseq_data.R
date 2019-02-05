setwd("D:/school_grad school/Project_borealis_sex_determination/Analysis/sample_clustering_analysis/PCA/")
#setwd("~/Desktop/borealis_sex_determination/PCA")

library(tidyverse)
library(graphics)

#another method Caroline recommended 
#dudi.pca
#library(mass)
#library(ad4)

#list of input_files
#input_file = "tropicalis_gonad.gene.counts.matrix"	
input_file = "tropicalis_gonad.gene.TMM.EXPR.matrix"
#input_file = "tropicalis_gonad.isoform.counts.matrix"
#input_file = "tropicalis_gonad.isoform.TMM.EXPR.matrix"

#read in the data
input_data <- read.table(input_file,header = TRUE, row.names ="Transcript_name")

#do some filtering
# I filtered out the genes with zeros expression in all samples;
#The second line is just there in case I more filtering to do
express_data_before <- input_data[rowSums(input_data)>1,]
express_data_head <- express_data_before[,]

#prcomp would expect the column to be variable (different genes) and row to be samples. Below is an example.
#Sample Gene1 Gene2 Gene3 Gene4
#Sample1  5 5 10`10`
#Sample2  2 3 15  29
#Hence I transcponse the data into the format that prcomp expects
#express_data <- log10(express_data)
express_data <- t(express_data_head)
express_data <- express_data[,colSums(express_data)>0]
express_data <- scale(express_data)

#Do PCA analysis using prcomp 
#By default, center = TRUE and it centers the variable to have mean equals to zero. 
#With parameter scale. = TRUE, we normalize the variables to have standard deviation equals to 1.
#prcomp return the following items:
###sdev: the standard deviations of the principal components
###rotation: the matrix of variable loadings (i.e., a matrix whose columns contain the eigenvectors). 
###x: if retx is true the value of the rotated data (the centred (and scaled if requested) data multiplied by the rotation matrix) is returned
###center, scale: the centering and scaling used
project_pca <- prcomp(express_data, retx=TRUE, center = TRUE, scale. = FALSE)

#calculate the project_pca_proportionvarianes
#way 1
project_pca_proportionvariances <- ((project_pca$sdev^2) / (sum(project_pca$sdev^2)))*100
#way 2
pca.var <- project_pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

#look at the pca result
#we would expect one PC for each gene; with this logic, we would the number of PC = number of genes
#but the number of samples (14 samples vs 246422 genes )put an upper bound on to the number of PC
#for detail explanation on why that is the case, please go to https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
summary(project_pca)

################################Visualization############################################
#visualize the pca rsult using pc.biplot -> this takes a lot of time to run, try to avoid using it if possible
#biplot(project_pca, choices = 2:3, scale = 1, pc.biplot = FALSE)

## plot pc1 and pc2
plot(project_pca$x[,1], project_pca$x[,2])

#Scree plot
#Scree plot allow you to look at the the percent of variation each PCA can account for
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")
barplot(project_pca_proportionvariances, cex.names=1, xlab=paste("Principal component (PC), 1-", length(project_pca$sdev)), ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))

## now make a fancy looking plot that shows the PCs and the variation:
library(ggplot2)

pca_plot <- data.frame(Sample=rownames(project_pca$x),
                       X=project_pca$x[,1],
                       Y=project_pca$x[,2])
pca_plot

ggplot(data=pca_plot, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", project_pca_proportionvariances[1], "%", sep="")) +
  ylab(paste("PC2 - ", project_pca_proportionvariances[2], "%", sep="")) +
  theme_classic() +
  ggtitle("X.tropicalis gonad RNAseq PCA") +
  geom_hline(mapping = NULL, data = NULL, yintercept = 0,
             na.rm = FALSE, show.legend = NA) +
  geom_vline(mapping = NULL, data = NULL, xintercept = 0,
             na.rm = FALSE, show.legend = NA)


#Pairs plots
par(cex=1.0, cex.axis=0.8, cex.main=0.8)
pairs(project_pca$x[,1:5], col="black", main="Principal components analysis bi-plot\nPCs 1-5", pch=16)

#Bi-plots
par(mar=c(4,4,4,4), mfrow=c(1,3), cex=1.0, cex.main=0.8, cex.axis=0.8)
#Plots scatter plot for PC 1 and 2
plot(project_pca$x, type="n", main="Principal components analysis bi-plot", xlab=paste("PC1, ", round(project_pca_proportionvariances[1], 2), "%"), ylab=paste("PC2, ", round(project_pca_proportionvariances[2], 2), "%"))
points(project_pca$x, col="blue", pch=16, cex=1)
#Plots scatter plot for PC 1 and 3
plot(project_pca$x[,1], project_pca$x[,3], type="n", main="Principal components analysis bi-plot", xlab=paste("PC1, ", round(project_pca_proportionvariances[1], 2), "%"), ylab=paste("PC3, ", round(project_pca_proportionvariances[3], 2), "%"))
points(project_pca$x[,1], project_pca$x[,3], col="black", pch=16, cex=1)

#Tri-plot
#This is very confusing to see, recommend not using it 
require(scatterplot3d)
par(mar=c(4,4,4,4), cex=1.0, cex.main=0.8, cex.axis=0.8)
scatterplot3d(project_pca$x[,1:3], angle=45, main="", color="black", pch=17, xlab=paste("PC1, ", round(project_pca_proportionvariances[1], 2), "%"), ylab=paste("PC2, ", round(project_pca_proportionvariances[2], 2), "%"), zlab=paste("PC3, ", round(project_pca_proportionvariances[3], 2), "%"), grid=FALSE, box=FALSE)

## get the name of the top 10 measurements (genes) that contribute
## most to pc1.
loading_scores <- project_pca$rotation[,1]
gene_scores <- abs(loading_scores) ## get the magnitudes
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_10_genes <- names(gene_score_ranked[1:10])

top_10_genes ## show the names of the top 10 genes

project_pca$rotation[top_10_genes,1] ## show the scores (and +/- sign)


###########################Other method#############################################
sd <- project_pca$sdev
loadings <- project_pca$rotation
rownames(loadings) <- colnames(express_data) 
scores <- project_pca$x 

#
quartz(height=7, width=7)
plot(scores[,1], scores[,2], xlab="PCA 1", ylab="PCA 2",
       type="n", xlim=c(min(scores[,1:2]), max(scores[,1:2])),
       ylim=c(min(scores[,1:2]), max(scores[,1:2])))
arrows(0,0,loadings[,1]*100,loadings[,2]*1000, length=0.1,
         angle=20, col="red") 

text(loadings[,1]*10*1.2,loadings[,2]*10*1.2,
     rownames(loadings), col="red", cex=0.7) 





data_for_PCA <- express_data %>%
  rownames_to_column %>% 
  gather(var, value, -rowname) %>% 
  spread(rowname, value) 
dim (data_for_PCA)

## calculate MDS (matrix of dissimilarities)
mds <- cmdscale(dist(data_for_PCA), eig=TRUE)  
mds$eig
