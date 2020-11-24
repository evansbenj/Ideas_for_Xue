#!/usr/bin/env Rscript

#-------------------------------------getting arguments from command line-------------------------------------
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = pasete(args[1], ".collapsed.tsv")
}


#-------------------------------------load package-------------------------------------
library(dplyr)
library(readr)
library(tidyr)
library(IRanges)
library(GenomicRanges)
#-------------------------------------extract gene coordinate from gff file-------------------------------------

# define input and output files
#gff_input = "project_sex_determination_gene/data/tropicalis_tad/genome/XENTR_9.1.gff3"
gff_input = "project_borealis_sexual_antagonism/data/genome/tropicalis/Xentr10.0.gene.gff3"

#read in mapping data
gff <- (read_tsv(gff_input, col_names = FALSE, skip =8) #read in data
        %>% filter(X3 == "gene")
        %>% select(chromosome =X1,
                   #feature = X3,
                   start = X4,
                   end = X5)
        %>% mutate (gene_position = paste(chromosome, ":", start, "-", end, sep = ""))
        # %>% extract(col=attribute, into = c("xenbase_ID", "Alias", "xenbase_name"),
        #             regex = "ID=([[:alnum:]]+);Alias=XT-9_1-([[:alnum:]]+);Name=([[:alnum:]]+)",
        #             remove = TRUE) #XENTR_9.1_gene_list.gff3
        # %>% select(-Alias)

)

#create GRange object for gff with ranges being genomic coordinate for each gene
gff_ir <- IRanges(start = gff$start, end = gff$end, names = gff$chromosome)
gff_gr <- GRanges(seqname=gff$chromosome,
                    ranges=gff_ir,
                    strand=NULL)

#-------------------------------------transcriptome mapping data-------------------------------------
# define input and output files
#mapping_input = "project_sex_determination_gene/data/tropicalis_tad/transcriptome/mapping_tropDNTrans_tropGenomeV91_gmap/tropDNTrans_tropGenomeV91_gmap_keep1match.tsv"
mapping_input = "project_sex_determination_gene/data/tropicalis_tad/mapping_tropDNTtrans_tropGenomeV10_gmap/tropDNTtrans_tropGenomeV10_gmap_keep1match.tsv"

#read in mapping data
mapping <- (read_tsv(mapping_input, col_names = FALSE) #read in data
            %>% select(trans_id = X1,
                       chromosome = X2,  
                       start = X3, 
                       end = X4) # name the columns 
            %>% mutate_at(vars(matches("chromosome")),list(factor))
            %>% mutate(width = end-start+1)
)

#create GRange object for chromosomes with ranges being genomic coordinate for each transcript
trans_ir <- IRanges(start = mapping$start, end = mapping$end, names = mapping$chromosome)
trans_gr <- GRanges(seqname=mapping$chromosome,
                       ranges=trans_ir,
                       strand=NULL)

#-------------------------------------gene-transcript overlap-------------------------------------
overlap <- findOverlaps(gff_gr, gff_gr, type = "any", ignore.strand=TRUE, minoverlap=1L, select = "all")

#extract matches 
gene_overlap_list <- ( tibble(gene_gene_number = queryHits(overlap), xenbase_gene_number =subjectHits(overlap))
                       %>% mutate(paper_gene_position = names(gene_paper_gr)[paper_gene_number],
                                  xenbase_gene_position = names(gene_xenbase_gr)[xenbase_gene_number])
                       %>% select(paper_gene_position, xenbase_gene_position)
)

