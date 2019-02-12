library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

#input_file_path <- file.path("D:", "school_grad school", 
input_file_path <- file.path("/Users/evanslab/Desktop",
                             "Xue_master_projects","Project_borealis_sex_determination/Analysis",
                             "sample_clustering_analysis")
output_file_path <- input_file_path
input_file <- paste(input_file_path, 
                    "identify_by_borealisGonad_DEgene",
                    "borealisGonadDEgene_4xBothSex_tropTran_blastnOut.tsv", 
                    sep="/")

input_data <- read_tsv(input_file, col_names = FALSE)

#select only the top trop hit for each borealise transcript
ortho_list <- input_data%>%
  group_by(X1)%>%
  arrange(X11, .by_group = TRUE)%>%
  top_n(1,X11)%>% #select the one with the lowest pvalue
  top_n(1,X12)%>% #if have ties for pvalue, select by bitscore
  ungroup() %>% 
  select(X2) %>% 
  rename("Transcript_name"="X2")%>%
  unique()

rm(input_data)

#read in the expression file
expression_file = paste(input_file_path,
                        "PCA/tropicalis_gonad.isoform.TMM.EXPR.matrix",
                        sep="/")
expression_data <- read_tsv(expression_file)

#select transcript that mapped to the borealise sex-biased gene
ortho_expression <- expression_data %>%
  inner_join(ortho_list, by=("Transcript_name"))%>%
  column_to_rownames("Transcript_name")
ortho_expression <- ortho_expression[rowSums(ortho_expression)>0,]
ortho_expression <- t(ortho_expression)
ortho_expression <- scale(ortho_expression)

####### Check scaling
#The below indicated that scaling did scale the data to have a mean of 0 and sd of 1
scaling_check <- as.data.frame (t(ortho_expression)) %>%
  rownames_to_column("Transcript_name")%>%
  gather(key = "samples", value = "expression_value", -Transcript_name)%>%
  group_by(Transcript_name)%>%
  mutate (Mean = mean(expression_value), SD =sd(expression_value) )

rm(scaling_check)

###### PCA
project_pca <- prcomp(ortho_expression, retx=TRUE, center = FALSE, scale. = FALSE)

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

##### kmeans
set.seed(200)
expCluster <- kmeans(ortho_expression, 2, nstart = 200)
expCluster$cluster

fviz_cluster(expCluster, data = ortho_expression)

