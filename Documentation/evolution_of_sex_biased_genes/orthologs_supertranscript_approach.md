directory: '/home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach
'
### check supertranscript
split and map then check if the parts mapped to same chromosome
```
cd /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach
perl ~/script/split_seq.pl SuperDuper.fasta > SuperDuper_splitted.fasta 

```
mapping: gmap
```
gmap -D /home/xue/borealis_transcriptome/borealis_genome_v1/db_borealis_genome_v1_gmap -d borealis_genome_v1 -A -B 5 -t 8 -f samse /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_unigeneApproach/borealis_laevis_transcriptome_laevis_unigene_besthit_unigeneSeq.fa | samtools view -S -b > /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_unigeneApproach/borealis_laevis_transcriptome_laevis_unigene_besthit_unigeneSeq_gmap.bam

```
