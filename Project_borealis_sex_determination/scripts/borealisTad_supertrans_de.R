library(tidyverse)
library(graphics)
library(dplyr)
library(edgeR)

#-------------------------------------mapping information-------------------------------------

#list of input_files for borealis

#input_file_folder <- file.path("Project_borealis_sex_determination/data/xb_transcriptome/supertrans/")
#input_file <- paste(input_file_folder,"borealisTad_denovoT_laevisGenome92_gmap_bedfile.bed", sep="/")

input_file_folder <- file.path("Project_borealis_sex_determination/data/analysis/transcriptome/mapping_trans_laevisGenome92_gmap")
input_file <- paste(input_file_folder,"borealisTad_denovoT_laevisGenome92_gmap_bedfile.bed", sep="/")





output_file_folder <- input_file_folder

#read in the data
input_data <- read_tsv(input_file, col_names = FALSE)

#filter the mapping and select one mapping per transcripts
mapping <- (input_data 
            %>% select(chromosome = X1,  
                      start = X2, 
                       end = X3,
                       supertrans_id = X4,
                       mapping_quality = X5,
                       strandness = X6) # name the columns 
            %>% group_by(supertrans_id) 
            %>% top_n(n=1, wt = mapping_quality) # only select the one with the highest mapping quality
            %>% add_count(supertrans_id, name = "count")
            %>% mutate(chrCount = sum(!grepl("Scaffold", chromosome) ),
                       scaCount = sum(grepl("Scaffold", chromosome)))
)

SingleMatch <- (mapping 
             %>% filter(count == 1)
             %>% select(supertrans_id, chromosome, start, end)
)


mulMatch_sameChr <- (mapping
             %>% group_by(supertrans_id)
             %>% filter(count >1, count==chrCount, n_distinct(chromosome)==1)
             %>% top_n(n=1, wt=start)
             %>% select(supertrans_id, chromosome, start, end)
)


mulMatch_allSca <- (mapping
                    %>% filter(count >1, count == scaCount)
                    %>% group_by(supertrans_id)
                    %>% top_n(n=1, wt=chromosome)
                    %>% select(supertrans_id, chromosome, start, end)
)

mulMatch_OneChr <- (mapping
                    %>% filter(count >1, chrCount == 1, count > chrCount, 
                               !grepl("Scaffold", chromosome))
                    %>% select(supertrans_id, chromosome, start, end)
)

filtered_mapping <- bind_rows(SingleMatch,mulMatch_sameChr, mulMatch_allSca, mulMatch_OneChr)

#-------------------------------------de--------------------------------------------
#read in the data

#supertrans
#input_file <- paste(input_file_folder,"borealis_tad_gonad_supertrans.counts.matrix", sep="/")

#transcriptome
input_file_folder <- file.path("Project_borealis_sex_determination/data/analysis/transcriptome/raw_count_kallisto")
input_file <- paste(input_file_folder,"borealis_tad_gonad_trans.counts.matrix", sep="/")

express_data <- read.table(input_file, header=T, row.names=1, com='')
express_data = express_data[rowSums(express_data)> 8,]


#stage information
stage46 <- c("XBO8","XBO15","XBO24","XBO27", "XBO12","XBO16","XBO23","XBO26")
stage48 <- c("XBO19","XBO20","XBO29","XBO30","XBO33","XBO34","XBO35","XBO36","XBO17","XBO28","XBO31","XBO32")
stage48_49 <- c("XBO19","XBO20","XBO21","XBO29","XBO30","XBO33","XBO34","XBO35","XBO36","XBO17","XBO28","XBO31","XBO32")

#sex information
stage46_sex <- factor(c(rep("female", 4), rep("male", 4)))
stage48_49_sex <- factor(c(rep("female", 8), rep("male", 5)))
stage48_sex <- factor(c(rep("female", 7), rep("male", 5)))

#select stage and conditions
rnaseqMatrix <- express_data[ ,which (names(express_data) %in% stage46)]
conditions = stage46_sex

# do DE with edgeR
exp_study = DGEList(counts=rnaseqMatrix, group=conditions)
exp_study = calcNormFactors(exp_study)
exp_study = estimateCommonDisp(exp_study)
exp_study = estimateTagwiseDisp(exp_study)
et = exactTest(exp_study, pair=c("female", "male"))
tTags = topTags(et,n=NULL)
result_table = tTags$table
result_table = data.frame(sampleA="female", sampleB="male", result_table)

#-------------------------------------combine DE and Mapping-------------------------------------

summary (result_table)

result_table <- rownames_to_column(result_table, var = "supertrans_id")
# result_filter <- (result_table
#                   %>% filter(FDR <= 0.05)
# 
# ) 

joined <- full_join(filtered_mapping,result_table,  by = "supertrans_id")

#-------------------------------------Summarize DE result by categories--------------------------
#extract chromosome and coordinate information 
FDR_cutoff = 0.05

de_table <- (joined 
  %>% filter (!is.na(start))
  %>% filter (!is.na(logFC))
  %>% mutate(category_broad = case_when (
        FDR > FDR_cutoff ~ "Non DE",
        chromosome == "chr8L" & start < 56690925 & end > 4605306 ~ "chr8L sex-link region",
        chromosome == "chr8L" ~ "chr8L non-sex-link region",
        chromosome == "chr8S" & start < 56690925 & end > 4605306 ~ "chr8S sex-link region",
        chromosome == "chr8S" ~ "chr8S non-sex-link region",
        grepl("Scaffold", chromosome) ~ "Scaffolds",
        grepl("MT", chromosome) ~ "Scaffolds",
        grepl("chr", chromosome) ~ "Other chromosomes"
      )
  )
  %>%mutate(category_broad_num = case_when(
        category_broad == "chr8L sex-link region" ~ "1",
        category_broad == "chr8L non-sex-link region" ~ "2",
        category_broad == "chr8S sex-link region" ~ "3",
        category_broad == "chr8S non-sex-link region" ~"4",
        category_broad == "Other chromosomes" ~ "5",
        category_broad == "Scaffolds" ~ "6",
        category_broad == "Non DE" ~ "7"
    )
  )
  %>% mutate(category_broad_cols = case_when (
        category_broad_num == "7" ~ "grey",
        category_broad_num == "6" ~ "blak",
        category_broad_num == "5" ~ "brown",
        category_broad_num == "4" ~ "deepskyblue1",
        category_broad_num == "3" ~ "blue",
        category_broad_num == "2" ~ "hotpink1",
        category_broad_num == "1" ~ "red1"
    )
  )
  %>% mutate(category_byChr = case_when (
        chromosome == "chr8L" & start < 56690925 & end > 4605306 ~ "chr8L-SL",
        chromosome == "chr8L" ~ "chr8L-NSL",
        chromosome == "chr8S" & start < 56690925 & end > 4605306 ~ "chr8S-SL",
        chromosome == "chr8S" ~ "chr8S-NSL",
        grepl("Scaffold", chromosome) ~ "Scaffolds",
        grepl("MT", chromosome) ~ "Scaffolds",
        TRUE ~ chromosome
    )
  )
)

de_summary_per_chr <- (de_table 
  %>%filter(!is.na(logFC))  #gene with very low raw count were not included in the DE analysis
  %>%mutate(
    DE_status = case_when(
      FDR < FDR_cutoff ~ "DE",
      FDR > FDR_cutoff ~ "Non-DE",
      )
    )
  %>%group_by(category_byChr)
  %>%count(category_byChr, DE_status)
  %>%rename(DE_count=n)
  %>%mutate(total_gene = sum(DE_count), proportion=DE_count/total_gene)
  %>%filter(DE_status=="DE")
)
#-------------------------------------binomial confidence interval-----------------------
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

labels = c("chr8L sex-link region", "chr8L non-sex-link region ", "chr8S sex-link region ", "chr8S non-sex-link region", "Other chromosomes", "Scaffolds (7.6%)", "Non DE")

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

