library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

input_file_path <- file.path("D:", "school_grad school", 
                             "Xue_master_projects/Project_borealis_sex_determination/Analysis",
                             "sample_clustering_analysis", "identify_by_borealisGonad_DEgene")

output_file_path <- input_file_path

input_file <- paste(input_file_path,"borealisGonadDEgene_tropTrna_blastnOut.tsv", sep="/")

input_data <- read_tsv(input_file, col_names = FALSE)

ortho_list <- input_data%>%
  group_by(X1)%>%
  arrange(X11, .by_group = TRUE)%>%
  top_n(1,X11)%>%
  top_n(1,X12)

rm(input_data)

ortho_trans <- ortho_list %>% 
  ungroup() %>% 
  select(X2) %>% 
  rename("Transcript_name"="X2")%>%
  unique()

expression_file = "../PCA/tropicalis_gonad.isoform.TMM.EXPR.matrix"

expression_data <- read.table(expression_file,header = TRUE, row.names ="Transcript_name")

ortho_expression <- expression_data %>%
  rownames_to_column("Transcript_name")%>%
  inner_join(ortho_trans, by=("Transcript_name"))%>%
  column_to_rownames("Transcript_name")

ortho_expression <- ortho_expression[rowSums(ortho_expression)>1,]
ortho_expression <- ortho_expression[,]
ortho_expression <- t(ortho_expression)
ortho_expression <- ortho_expression[,colSums(ortho_expression)>0]
ortho_expression <- scale(ortho_expression)
project_pca <- prcomp(ortho_expression, retx=TRUE, center = TRUE, scale. = FALSE)

#calculate the project_pca_proportionvarianes
#way 1
project_pca_proportionvariances <- ((project_pca$sdev^2) / (sum(project_pca$sdev^2)))*100
#way 2
pca.var <- project_pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

## now make a fancy looking plot that shows the PCs and the variation:

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

#kmeans
set.seed(20)
expCluster <- kmeans(ortho_expression, 2, nstart = 20)
expCluster$cluster

fviz_cluster(expCluster, data = ortho_expression)

