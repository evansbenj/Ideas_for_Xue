# Identification of orthologs by mapping all the transcripts to laevis genome
## existing mapping information
```
```
## more mapping
mapping tropicalis transcriptome to laevis genome
```
#gmap: mapping tropicalis **de novo** transcriptome to laevis genome
gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -A -Z -f samse /home/xue/borealis_DE/liver_mvsf/filtered_edgeRout/liver_trans_fdr005_header.fa > /home/xue/borealis_DE/liver_mvsf/mapping_GMAP/liver_DE_gmap_out.sam


#gmap: mapping tropicalis genome guided transcriptome to laevis genome
gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -A -Z -f samse /home/xue/borealis_DE/liver_mvsf/filtered_edgeRout/liver_trans_fdr005_header.fa > /home/xue/borealis_DE/liver_mvsf/mapping_GMAP/liver_DE_gmap_out.sam


```
