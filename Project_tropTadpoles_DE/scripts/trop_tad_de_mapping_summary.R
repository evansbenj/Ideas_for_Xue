library(tidyverse)
library(edgeR)
library(IRanges)
library(GenomicRanges)

#set some cutoff
FDR_cutoff = 0.1
logFC_cutoff_up = 1
logFC_cutoff_down = -1

#-------------------------------------read in transcriptome de data-------------------------------------
#read in de data
input_file_folder <- file.path("project_sex_determination_gene/data/tropicalis_tad/transcriptome/trop_tad_de_result")
input_file <- paste(input_file_folder,"de_result_trop_stage50.tsv", sep="/")

de_data <- read.table(input_file, header=T, row.names = 1, com='')

de_data <- rownames_to_column(de_data, var = "transcript_id")

#-------------------------------------gene list from paper-------------------------------------
input_file_folder <- file.path("project_sex_determination_gene/data/tropicalis_tad/transcriptome/gene_sexlink_region")
input_file <- paste(input_file_folder,"sex_gene_list_from_paper.csv", sep="/")

gene_list <- read_csv(input_file)

gene_position <- (gene_list
                  %>% select(gene_position = "gene position")
                  %>% extract(col=gene_position, into = c("chromosome", "start", "end"), 
                              regex = "([[:alnum:]]+):([[:digit:]]+)-([[:digit:]]+)", 
                              remove = FALSE)
                  %>% filter(!is.na(start))
                  %>% mutate(chromosome = as.factor(chromosome),
                             start = as.numeric(start), 
                             end = as.numeric(end))
                  #%>% select(chromosome, start, end, gene_position, )
)


#-----------------------------------Genome mapping - identify gene and corresponding transcripts by overlapping ranges-------------------------------------
#read in genome mapping data
input_file_folder <- file.path("project_sex_determination_gene/data/tropicalis_tad/transcriptome/mapping_tropDNTrans_tropGenomeV91_gmap")
input_file <- paste(input_file_folder,"tropDNTrans_tropGenomeV91_gmap_keep1match.tsv", sep="/")

genome_mapping_data <- (read_tsv(input_file, col_names = FALSE)
                 %>% select(transcript_id = X1,
                            chromosome = X2,
                            start = X3,
                            end = X4)
                 %>% mutate(chromosome = as.factor(chromosome))
                 %>% mutate(genome_alignment_position = paste(chromosome, ":", start, "-", end, sep = ""))
)

#just want to look at how many transcript mapped to chr07
trans_list_chr07 <- (left_join(de_data, genome_mapping_data, by= "transcript_id")
                     %>% filter(chromosome == "Chr07")
                     #%>% filter(logFC>=logFC_cutoff_up | logFC <= logFC_cutoff_down)
                     %>% filter(FDR <=FDR_cutoff)
)



#create GRange object for GENE list
gene_ir <- IRanges(start = gene_position$start, end = gene_position$end, names = gene_position$gene_position)
gene_gr <- GRanges(seqname=gene_position$chromosome,
                       ranges=gene_ir,
                       strand=NULL)

#create GRange object for TRANSCRIPT list
trans_ir <- IRanges(start = trans_list_chr07$start, end = trans_list_chr07$end, names = trans_list_chr07$transcript_id)
trans_gr <- GRanges(seqname = trans_list_chr07$chromosome,
                       ranges=trans_ir,
                       strand=NULL)

#match gene list with transcripts list
overlap <- findOverlaps(gene_gr, trans_gr, type = "any", ignore.strand=TRUE, minoverlap=1L, select = "all")

#extract matches 
gene_transcript_list <- ( tibble(gene_number = queryHits(overlap), transcript_number =subjectHits(overlap))
                         %>% mutate(gene_position = names(gene_gr)[gene_number],
                                    transcript_id = names(trans_gr)[transcript_number])
                         %>%select(gene_position, transcript_id)
                         %>%left_join(trans_list_chr07, by =c("transcript_id"="transcript_id"))
                         %>%select(gene_position, transcript_id, logFC, FDR, genome_alignment_position)
                         %>% filter(gene_position == "Chr07:663201-668367" | gene_position =="Chr07:1928366-1929352"|
                                      gene_position == "Chr07:1376005-1401810"|
                                      gene_position == "Chr07:1927565-1928168"|
                                      gene_position == "Chr07:1936551-1940958"| gene_position =="Chr07:2045094-2046084"|
                                      gene_position == "Chr07:3299492-3300392")
                         #%>% group_by(gene_position)
)

output_folder <- file.path("project_sex_determination_gene/data/tropicalis_tad/transcriptome/gene_sexlink_region")
output_file <- paste(output_folder,"gene_de_transcript_list.tsv", sep="/")
write_tsv(gene_transcript_list, output_file, na = "NA", append = FALSE, col_names = FALSE,
          quote_escape = "none")


#output gene_position in a bed file format
output_file <- paste(output_folder,"gene_position.bed", sep="/")
write_tsv(gene_position, output_file, na = "NA", append = FALSE, col_names = FALSE,
          quote_escape = "none")

#-----------------------------------------Gene mapping - mapping transcripts to the gene sequences----------------------
#read in mapping data - transcriptome mapped to the 50 gene sequences
input_file_folder <- file.path("project_sex_determination_gene/data/tropicalis_tad/transcriptome/gene_sexlink_region")
input_file <- paste(input_file_folder,"tropTadGonad_DNTrans_gene_fromPaper_gmap_keep1match.tsv", sep="/")

gene_mapping_data <- (read_tsv(input_file, col_names = FALSE)
                 %>% select(transcript_id = X1,
                            gene = X2,
                            start = X3,
                            end = X4)
                 %>% extract(col=gene, into = c("gene_chromosome", "gene_start", "gene_end"), 
                                 regex = "([[:alnum:]]+)_([[:digit:]]+)-([[:digit:]]+)", 
                                 remove = FALSE)
                 %>% mutate(gene = paste0(gene_chromosome,":",gene_start,"-",gene_end),
                            alignment_length = end-start+1, 
                            gene_length = as.numeric(gene_end) - as.numeric(gene_start) +1)
                 %>% select(-gene_chromosome, -gene_start, -gene_end)
                 #%>% filter(alignment_length > 200)
                 #%>% select(transcript_id, gene_position)
                 %>% distinct()
                 %>% left_join(de_data, by = "transcript_id")
                 # %>% filter(gene == "Chr07:663201-668367" | gene =="Chr07:1928366-1929352"|
                 #              gene == "Chr07:1376005-1401810"| 
                 #              gene == "Chr07:1927565-1928168"| 
                 #              gene == "Chr07:1936551-1940958"| gene =="Chr07:2045094-2046084"|
                 #              gene == "Chr07:3299492-3300392")
                 %>% filter(!is.na(FDR), logFC > logFC_cutoff_up|logFC < logFC_cutoff_down, FDR < 0.5)
                 #%>% group_by(gene)
                 #%>% top_n(5, wt = alignment_length)
                 %>% left_join(mapping_data, by = "transcript_id")
                 %>% filter(grepl (pattern = "caf", chromosome)|chromosome == "Chr07")
                 %>%left_join(gene_list, by =c(gene = "gene position"))
                 %>% filter (gene_position!="Chr07:1363009-1369803"| gene_position!="Chr07:1685765-1691609" )
                 

)


#-----------------------------------------compile and output table----------------------
combine_list <- (left_join(gene_transcript_list, gene_mapping_data, by = "transcript_id" )
                 %>% filter (gene_position != gene)
                 %>%left_join(gene_list, by =c(gene_position = "gene position"))
                 %>%select(gene_name = "gene name", gene_position, gene_length,
                           transcript_id, logFC, FDR,
                           GeneMapping_align_length=alignment_length)
)


output_folder <- file.path("project_sex_determination_gene/data/tropicalis_tad/transcriptome/gene_sexlink_region")
output_file <- paste(output_folder,"gene_transcript_list.tsv", sep="/")
write.table(combine_list, output_file, sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE,)





