library("limma")
library("edgeR")
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyverse)

species <- "borealis"
tissue <- "liver"
quantification <-"mappedOnly"
rawcount_cutoff <- 8
FDR_cutoff <- 0.05


num_female_samples <- 4
num_male_sample <- 4

#select input file based on specis and tissue

input_file <- "Project_borealis_sexual_antagonism/data/supertranscriptome/counts_salmon/borealis_liver_supertrans_counts.matrix"

data = read.table(input_file, header=T, row.names=1, com='')
col_ordering = c(1,2,3,4,5,6,7,8)
rnaseqMatrix = data[,col_ordering]

samps <- read.csv(file.path("Project_borealis_sexual_antagonism/data/supertranscriptome/DTU_analysis", "sample_condition.csv"))
samps$condition <- factor(samps$condition) 
table(samps$condition)
conditions = factor(c(rep("female_rep", 4), rep("male_rep", 4)))


#different filtering method
#filtering method 1
rnaseqMatrix = rnaseqMatrix[rowSums(rnaseqMatrix)>=rawcount_cutoff,] 
#filtering method 2
#rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 0) >= 1,]
#filtering method 3
# cpm_log <- cpm(rnaseqMatrix, log = TRUE)
# median_log2_cpm <- apply(cpm_log, 1, median)
# expr_cutoff <- -3
# rnaseqMatrix <- rnaseqMatrix[median_log2_cpm > expr_cutoff, ]

exp_study = DGEList(counts=rnaseqMatrix, group=conditions)
exp_study = calcNormFactors(exp_study)
exp_study = estimateCommonDisp(exp_study)
exp_study = estimateTagwiseDisp(exp_study)
et = exactTest(exp_study, pair=c("female_rep", "male_rep"))
tTags = topTags(et,n=NULL)
result_table = tTags$table
result_table = data.frame(sampleA="female", sampleB="male", result_table)


#-----------------------------------extract coordinat information ---------------------
bed_file <- "Project_borealis_sexual_antagonism/data/supertranscriptome/mapping_laevis_genomeV92_gamp/borealis_superTrans_laevisV92_genome_gmap_bedfile.bed"

bed_data <- read_tsv(bed_file, col_names =FALSE)

bed_filtered <- (bed_data 
                 %>% group_by(X4)
                 %>% arrange(X4,desc(X5))
                 %>% top_n(1,X5)
                 %>% rename(Chr = X1, Start =X2, End=X3,transcript_id=X4)
                 %>% select(-X5, -X6)
                 #%>% column_to_rownames(var="transcript_id")
                 )

result_table <- result_table%>% rownames_to_column(var = "transcript_id")
de_table <- inner_join(result_table, bed_filtered, by="transcript_id")

#---------------------------------------Summarize DE result by categories--------------------------
#extract chromosome and coordinate information 
#turn the data.frame into tibble and categorize the genes based on genomic locations
de_table <- de_table %>%
  mutate( 
    category_broad = case_when (
      FDR > FDR_cutoff ~ "Non DE",
      Chr == "chr8L" & Start < 56690925 & End > 4605306 ~ "chr8L sex-link region",
      Chr == "chr8L" ~ "chr8L non-sex-link region",
      Chr == "chr8S" & Start < 56690925 & End > 4605306 ~ "chr8S sex-link region",
      Chr == "chr8S" ~ "chr8S non-sex-link region",
      grepl("Scaffold", Chr) ~ "Scaffolds",
      grepl("MT", Chr) ~ "Scaffolds",
      grepl("chr", Chr) ~ "Other chromosomes"
    )
  )%>%
  mutate(
    category_broad_num = case_when(
      category_broad == "chr8L sex-link region" ~ "1",
      category_broad == "chr8L non-sex-link region" ~ "2",
      category_broad == "chr8S sex-link region" ~ "3",
      category_broad == "chr8S non-sex-link region" ~"4",
      category_broad == "Other chromosomes" ~ "5",
      category_broad == "Scaffolds" ~ "6",
      category_broad == "Non DE" ~ "7"
    )
  )%>%
  mutate(
    category_broad_cols = case_when (
      category_broad_num == "7" ~ "grey",
      category_broad_num == "6" ~ "blak",
      category_broad_num == "5" ~ "brown",
      category_broad_num == "4" ~ "deepskyblue1",
      category_broad_num == "3" ~ "blue",
      category_broad_num == "2" ~ "hotpink1",
      category_broad_num == "1" ~ "red1"
    )
  )%>%
  mutate( 
    category_byChr = case_when (
      Chr == "chr8L" & Start < 56690925 & End > 4605306 ~ "chr8L-SL",
      Chr == "chr8L" ~ "chr8L-NSL",
      Chr == "chr8S" & Start < 56690925 & End > 4605306 ~ "chr8S-SL",
      Chr == "chr8S" ~ "chr8S-NSL",
      grepl("Scaffold", Chr) ~ "Scaffolds",
      grepl("MT", Chr) ~ "Scaffolds",
      TRUE ~ Chr
    )
  )

total_DE_count <- sum(de_table$FDR < FDR_cutoff) 
de_summary_category_broad <- de_table %>%filter(category_broad!="Non DE")%>%group_by(category_broad) %>%summarise (count = n(), percent = n()/total_DE_count*100)

de_summary_per_chr_gene <- data.frame (de_table %>% group_by(category_byChr) %>%
                                    summarize (total_gene =n())) 

de_summary_per_chr_de <- de_table %>%
  mutate(
    category_byChr = case_when (
      FDR > 0.05 ~ "non-DE",
      TRUE ~ category_byChr
    )  
  ) %>%
  group_by(category_byChr) %>%
  summarise (DE_count = n())

de_summary_per_chr <- merge.data.frame(de_summary_per_chr_gene, de_summary_per_chr_de, by = c("category_byChr"), all = TRUE)

de_summary_per_chr <- de_summary_per_chr[-c(21), ]
de_summary_per_chr[is.na(de_summary_per_chr)] <- 0

#-------------------------------------------binomial confidence interval-----------------------
de_summary_per_chr_ci <- de_summary_per_chr %>% 
                          mutate(proportion=DE_count/total_gene) %>% 
                          mutate(binomial_variance = total_gene*proportion*(1-proportion))%>%
                          mutate(upper=DE_count + sqrt(binomial_variance)*1.96) %>%
                          mutate(lower=DE_count - sqrt(binomial_variance)*1.96) %>%
                          mutate(uppercentile = upper/total_gene)%>%
                          mutate(lowercentile = lower/total_gene)%>%
                          mutate( 
                            chr_cols = case_when (
                              category_byChr == "chr8L-SL" ~ "red1",
                              category_byChr == "chr8L-NSL" ~ "hotpink1",
                              category_byChr == "chr8S-SL" ~ "blue",
                              category_byChr == "chr8S-NSL" ~"deepskyblue1",
                              category_byChr == "chr2L" ~ "purple",
                              TRUE ~ "black"
                            )
                          )

col <- as.character(de_summary_per_chr_ci$chr_cols)
names(col) <- as.character(de_summary_per_chr_ci$category_byChr)

de_ci_plot <- ggplot(de_summary_per_chr_ci, aes(x=category_byChr, y=proportion, label = category_byChr))+
  geom_point(aes(size=3,colour = category_byChr), show.legend = FALSE)+
  geom_errorbar(aes(ymin=uppercentile, ymax=lowercentile), width=0.4,colour = col)+
  xlab("Chromosomes")+
  ylab("% sex-biased transcripts per total genes")+
  scale_x_discrete(drop = FALSE)+
  scale_colour_manual(values = col)+
  theme_classic()+
  theme(
    axis.title.y = element_text(size = 20, margin = margin(0,20,0,0)),
    axis.title.x = element_text(size = 20, margin = margin(0,20,0,0)),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", colour = col),
    axis.text = element_text(size = 18)
  ) 

de_ci_plot

#------------------------------------make a pretty MA plot for the DE result-------------------------------------
ma_col <- c("1" = "red1", "2" = "hotpink2", "3" = "blue", "4" = "steelblue3", "5" = "black", "6" = "brown4", "7" = "grey")

ma_legend_label <- 
labels = c("chr8L sex-link region (1.7%)", "chr8L non-sex-link region (6.8%)", "chr8S sex-link region (0.4%)", "chr8S non-sex-link region (6.4%)", "Other chromosomes (77.1%)", "Scaffolds (7.6%)", "Non DE")

#to make the point transparent (the transparency level is known as the alpha), the lower the alpha, the more transparent the point will be
de_table<- de_table %>% arrange(desc(category_broad_num))
ma_plot<-ggplot(de_table) + 
      geom_point(
        aes(x=logCPM, y =logFC,col = category_broad_num), 
        size=3, alpha=0.6, show.legend =TRUE
      ) +
      xlab("logCPM") + 
      ylab("log2FC") + 
      guides(alpha=FALSE, size = FALSE)+ #don't include alpha and size infor on legend
      scale_colour_manual( #manipulate the color of the points and legends
        values = ma_col, #tell ggplot what color to use for the points 
        name = "Genomic location", #name of the legend
        drop = FALSE,
        labels = labels
      )+
      theme_classic() +
      theme(
        axis.title.y = element_text(size=20, margin=margin(0,20,0,0)),
        axis.title.x = element_text(size=20, margin=margin(20,0,0,0)),
        axis.text = element_text(size = 18),
        legend.title = element_text(colour="black", size=13, face="bold"),
        legend.text = element_text(colour="black", size = 13),
        legend.justification = "top",
        legend.position = c(0.85, 1),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")
      )
ma_plot
#-------------------------------------------------output data-----------------------
#output file names
output_prefix = paste(species,"_", tissue,"_",quantification,"_rawcount",rawcount_cutoff, "_FDR",(FDR_cutoff*100), sep="")

#output DE result - unfiltered
write.table(result_table, file=(paste(output_prefix,"_deFullResult.tsv", sep="")), sep='\t', quote=F, row.names=T)

#output summary (by chr) of DE result 
write.table(de_summary_per_chr_ci, file=(paste(output_prefix,"_SummaryPerChr.csv", sep=",")), sep=',', quote=F, row.names=F)

#output the transcript ID of DE in the 8L sex-link region
write.table(de_table[which(de_table$category == "chr8L sex-link region"),], file=(paste(output_prefix,"_8LSexRegion_TransID.csv", sep="")),sep=',', quote=F, row.names=T)

#distribution plot with binomial confidence interval
ggsave(filename=(paste(output_prefix,"_distributionPlot_withBionomialCI.pdf", sep="")), plot=de_ci_plot, height = 8.27, width = 11.69,units = c("in"), device = "pdf")

#transcript bin expression
ggsave(filename=(paste(output_prefix,"_MA_plot.pdf", sep="")), plot=ma_plot, height = 8.27, width = 11.69,units = c("in"), device = "pdf")

dev.off()