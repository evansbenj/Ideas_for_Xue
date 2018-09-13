mapping super transcript to laveis genome
```
gmap -d laevis92_gmap/ -D /home/xue/genome_data/laevis_genome/db_gmap_xl92 -A -B 5 -t 24 -f samse --cross-species /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_supertranscriptApproach/SuperDuper.fasta | samtools view -S -b > /home/xue/borealis_transcriptome/borealis_supertranscriptome_2018/SuperDuper_gmap_laevisGenome.bam

```
