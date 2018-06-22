# Identification of orthologs by mapping all the transcripts to laevis genome
## existing mapping information
```
#borealis de transcript to laevis genome v92
/home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_byGenomicLocation/mapped_gmap/liver_DE_gmap_out_bedfile.tsv

#laevis denovo transcriptome to laevis genome v92
/home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_byGenomicLocation/mapped_gmap/laevis_denovo_transcriptome_genome_gmap_bedfile.tsv
```
## Mapping
I have done some of trancriptome-laevis genome mapping for laevis denovo transcriptome and borealis de transcript before but didn't do it for tropicalis transcriptome and borealis transcriptome. I dont need the mapping information for the borealis transcriptome-laeivs genome for now but will need it for future work.  
- mapping the entire borealis transcriptome to laevis genome
  ```
  #gmap: mapping borealis **de novo** transcriptome to laevis genome
  gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -A -B 5 -t 8 -f samse --cross-species /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta | samtools view -S -b > /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_laevisGenomeApproach/borealis_denovoT_laevisV92_genome_gmap.bam

  #gmap: mapping borealis genome guided transcriptome to laevis genome
  gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -A -B 5 -t 8 -f samse --cross-species /home/xue/borealis_transcriptome/borealis_gg_transcriptome_may2018/borealis_ggT_trinityOut.fasta | samtools view -S -b > /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_laevisGenomeApproach/borealis_ggT_laevisV92_genome_gmap.bam
  ```

  - mapping tropicalis transcriptome to laevis genome
  ```
  #gmap: mapping tropicalis **de novo** transcriptome to laevis genome
  gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -A -B 5 -t 8 -f samse --cross-species /home/xue/tropicalis_transcriptome/tropicalis_denovo_transcriptome/tropicalis_trinityout.Trinity.fasta | samtools view -S -b > /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/laevis_tropicalis_orthologs_laevisGenomeApproach/tropicalis_denovoT_laevisV92_genome_gmap.bam

  #gmap: mapping tropicalis genome guided transcriptome to laevis genome
  gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -A -B 5 -t 8 -f samse --cross-species /home/xue/tropicalis_transcriptome/tropicalis_gg_transcriptome/tropicalis_gg_transcriptome.fasta | samtools view -S -b > /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/laevis_tropicalis_orthologs_laevisGenomeApproach/tropicalis_ggT_laevisV92_genome_gmap.bam
  ```

## extract laevis gene information from laevis gff3 file
```
awk '$3 == "gene"{print}' XENLA_9.2_Xenbase.gff3 > XENLA_gene.gff3
```
