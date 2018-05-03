# path to files:
```
#path to borealis_genome_v1 from Austin


```

# Trinity: genome guided asssembly
```
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq  --left //home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/laevis_liver_R1_scythe.fastq.gz --right /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/laevis_liver_R2_scythe.fastq.gz --genome_guided_bam /home/xue/laevis_transcriptome_mar2018/laevis_RNAseq_star/laevis_RNAseq_StarAligned.sortedByCoord.out.bam --genome_guided_max_intron 20000 --CPU 25 --full_cleanup --max_memory 200G --output /home/xue/laevis_transcriptome_mar2018/laevis_transcriptome_trinityout; echo "trinity is done at info115 in screen laevis" | mail sarahsongxy@gmail.com

```
