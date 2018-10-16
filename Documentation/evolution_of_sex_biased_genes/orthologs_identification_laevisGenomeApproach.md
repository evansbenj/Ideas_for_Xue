# Identification of orthologs by mapping all the transcripts to laevis genome
## Mapping
I have done some of trancriptome-laevis genome mapping for laevis denovo transcriptome and borealis de transcript before but didn't do it for tropicalis transcriptome and borealis transcriptome. I dont need the mapping information for the borealis transcriptome-laeivs genome for now but will need it for future work.  
- mapping the entire borealis transcriptome to laevis genome
  ```
  #gmap: mapping borealis **de novo** transcriptome to laevis genome
  gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -A -B 5 -t 8 -f samse --cross-species /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta | samtools view -S -b > /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_laevisGenomeApproach/borealis_denovoT_laevisV92_genome_gmap.bam

  #gmap: mapping borealis genome guided transcriptome to laevis genome
  gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -A -B 5 -t 8 -f samse --cross-species /home/xue/borealis_transcriptome/borealis_gg_transcriptome_may2018/borealis_ggT_trinityOut.fasta | samtools view -S -b > /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_laevisGenomeApproach/borealis_ggT_laevisV92_genome_gmap.bam
  
  #the transcriptome are too large, so mapping the entire transcriptome to the genome fail after running gmap for two weeks. I am going to split the transcriptome into separate file, each contains 100000 transcripts
   cd /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017 
  time gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -A -B 5 -t 8 -f samse --cross-species /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/subset_15.fasta | samtools view -S -b > /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/borealis_denovoT_laevisV92_genome_gmap_subset15.bam
  ```
- mapping tropicalis transcriptome to laevis genome
  ```
  #gmap: mapping tropicalis **de novo** transcriptome to laevis genome
  gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -A -B 5 -t 8 -f samse --cross-species /home/xue/tropicalis_transcriptome/tropicalis_denovo_transcriptome/tropicalis_trinityout.Trinity.fasta | samtools view -S -b > /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/laevis_tropicalis_orthologs_laevisGenomeApproach/tropicalis_denovoT_laevisV92_genome_gmap.bam

  #gmap: mapping tropicalis genome guided transcriptome to laevis genome
  gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -A -B 5 -t 8 -f samse --cross-species /home/xue/tropicalis_transcriptome/tropicalis_gg_transcriptome/tropicalis_gg_transcriptome.fasta | samtools view -S -b > /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/laevis_tropicalis_orthologs_laevisGenomeApproach/tropicalis_ggT_laevisV92_genome_gmap.bam
  ```
## path to mapping information: all done by GMAP
```
#borealis de transcript to laevis genome v92
/home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_byGenomicLocation/mapped_gmap/liver_DE_gmap_out.bam

#laevis denovo transcriptome to laevis genome v92
/home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_byGenomicLocation/mapped_gmap/laevis_denovo_transcriptome_genome_gmap_bedfile.tsv

#tropicalis **de novo** transcriptome to laevis genome v92
/home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/laevis_tropicalis_orthologs_laevisGenomeApproach/tropicalis_denovoT_laevisV92_genome_gmap.bam
```

## filter bam file
filter the mapping read remove the unmapped reads. Then I pile the filtered output to bedtools, which extracts alignment coordinates for transcripts based on their CIGAR strings;
```
#working directory
cd /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_laevisGenomeApproach
#move everything in the laevis_tropicalis_orthologs_laevisGenomeApproach to here so that I can have all the file in one place
mv ../laevis_tropicalis_orthologs_laevisGenomeApproach/* .

#filtering mapping outputs and extract coordinate
samtools view -F 0x04 -b /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_byGenomicLocation/mapped_gmap/liver_DE_gmap_out.bam | bedtools bamtobed -i > /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_laevisGenomeApproach/borealis_de_laevisV92_genome_gmap_bedfile.tsv

samtools view -F 0x04 -b /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_byGenomicLocation/mapped_gmap/laevis_denovo_transcriptome_genome_gmap.bam | bedtools bamtobed -i >  /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_laevisGenomeApproach/laevis_denovoT_laevisV92_genome_gmap_bedfile.tsv

samtools view -F 0x04 -b tropicalis_denovoT_laevisV92_genome_gmap.bam | bedtools bamtobed -i > tropicalis_denovoT_laevisV92_genome_gmap_bedfile2.tsv
```
## filter the bedfile and keep only the top hit
filter the bedfile and keep only the hit with the highest mapping quality score; if there are two hit with the same mapping quality score, keep the first one only;
```
perl ~/script/orthologs_identification_filterBedfile.pl borealis_de_laevisV92_genome_gmap_bedfile.bed > borealis_de_laevisV92_genome_gmap_bedfile_filtered.tsv

perl ~/script/orthologs_identification_filterBedfile.pl laevis_denovoT_laevisV92_genome_gmap_bedfile.bed > laevis_denovoT_laevisV92_genome_gmap_bedfile_filtered.tsv

perl ~/script/orthologs_identification_filterBedfile.pl tropicalis_denovoT_laevisV92_genome_gmap_bedfile.bed > tropicalis_denovoT_laevisV92_genome_gmap_bedfile_filtered.tsv
```

## extract laevis gene information from laevis gff3 file
```
cd /home/xue/genome_data/laevis_genome/gff3_xl9_2/
awk '$3 == "gene"{print}' XENLA_9.2_Xenbase.gff3 > /home/xue/genome_data/laevis_genome/gff3_xl9_2/XENLA_gene.gff3

#extract the sequence of each gene

#blastn gff gene to laevis genome
time blastn -task megablast -db ~/genome_data/laevis_genome/db_blastn_laevisGenome/Xl9_2_blastn_db -outfmt 6 -evalue 0.00005 -query /home/xue/genome_data/laevis_genome/gff3_xl9_2/XENLA_gene.fa -out /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_laevisGenomeApproach/laevis_genes_laevisV92_genome_blastout.tsv

time blastn -task blastn -db ~/genome_data/laevis_genome/db_blastn_laevisGenome/Xl9_2_blastn_db -outfmt 6 -evalue 0.00005 -query /home/xue/genome_data/laevis_genome/gff3_xl9_2/XENLA_gene.fa -out /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_laevisGenomeApproach/laevis_genes_laevisV92_genome_blastnout.tsv

#blast expanded_part_borealis_mapped_uncollapsed_transcript gff genome to laevis genome
time blastn -task blastn -db ~/genome_data/laevis_genome/db_blastn_laevisGenome/Xl9_2_blastn_db -outfmt 6 -evalue 0.00005 -query /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/expanded_laevis_gff_genome/expanded_laevis_gffGenome_laevis_genome_blastn/subset_15.fasta -out /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/expanded_laevis_gff_genome/expanded_laevis_gffGenome_laevis_genome_blastn/subset_15_blastout
```



## check for overlap between laevis gff genes and transcripts
```
bedtools intersect -wao -a /home/xue/genome_data/laevis_genome/gff3_xl9_2/XENLA_gene.gff3 -b borealis_de_laevisV92_genome_gmap_bedfile.bed, laevis_denovoT_laevisV92_genome_gmap_bedfile.bed, tropicalis_denovoT_laevisV92_genome_gmap_bedfile.bed > orthologs_list.tsv
```
