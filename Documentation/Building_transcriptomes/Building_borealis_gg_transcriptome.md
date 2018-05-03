# path to files:
```
#path to borealis_genome_v1 from Austin
/home/xue/borealis_transcriptome/borealis_genome_v1/Xbo.v1.fa.gz
#path to the genome-guided transcriptome
```

# Trinity: genome-guided transcriptome assembly
```
#Pre-assembly: mapping reads to borealis genome and create a genome_guided_bam file using STAR
STAR --runThreadN 20 --genomeDir ~/genome_data/laevis_genome/db_star_laevisGenome_wGTF --readFilesCommand zcat --readFilesIn /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/laevis_liver_R1_scythe.fastq.gz /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/laevis_liver_R2_scythe.fastq.gz --outFileNamePrefix /home/xue/laevis_transcriptome_mar2018/laevis_RNAseq_Star --outSAMtype BAM SortedByCoordinate

#Trinty assembly with 
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq  --left //home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/laevis_liver_R1_scythe.fastq.gz --right /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/laevis_liver_R2_scythe.fastq.gz --genome_guided_bam /home/xue/laevis_transcriptome_mar2018/laevis_RNAseq_star/laevis_RNAseq_StarAligned.sortedByCoord.out.bam --genome_guided_max_intron 20000 --CPU 25 --full_cleanup --max_memory 200G --output home/xue/borealis_transcriptome/borealis_gg_transcriptome_may2018/borealis_gg_trinityout; echo "trinity is done at info115 in screen laevis" | mail sarahsongxy@gmail.com

```
