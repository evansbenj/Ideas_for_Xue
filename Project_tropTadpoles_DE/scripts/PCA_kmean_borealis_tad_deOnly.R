#pull out DE and do a pca

library(tidyverse)
library(graphics)


#input_files

#supertrans
#input_file = "Project_borealis_sex_determination/data/borealis_tad/supertrans/raw_count_kallisto/borealis_tad_gonad_supertrans.TMM.EXPR.matrix"
#trans
input_file = "Project_borealis_sex_determination/data/borealis_tad/transcriptome/raw_count_kallisto/borealis_tad_gonad_trans.TMM.EXPR.matrix"

#read in the data
input_data <-  read.table(input_file, header=T, row.names=1, com='')

#stage information
stage46 <- c("XBO8","XBO15","XBO24","XBO27", "XBO12","XBO16","XBO23","XBO26")
stage48_49 <- c("XBO19","XBO20","XBO21","XBO29","XBO30","XBO33","XBO34","XBO35","XBO36","XBO17","XBO28","XBO31","XBO32")
DE_list <- (de_table_46t 
            %>% filter(!is.na(logFC) && FDR > 0 && FDR <=0.05 )
            %>% select(supertrans_id)
)

#select data
express_data <- input_data[rowSums(input_data) > 1,] #sum of count > 1
express_data <- express_data[,which (names(express_data) %in% stage46)] 
express_data <- express_data[which (row.names(express_data) %in% DE_list$supertrans_id),] 

################################PCA############################################
#prcomp would expect the column to be variable (different genes) and row to be samples.
express_data <- t(express_data)
express_data <- express_data[,colSums(express_data)>0]

#Do PCA analysis using prcomp 
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
