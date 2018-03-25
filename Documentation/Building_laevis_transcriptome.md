

# Building laevis transcriptome with Trinity (genome-guided method)

### Pre-assembly: mapping reads to laevis genome
Time cost:  92min
```
STAR --runThreadN 20 --genomeDir ~/genome_data/laevis_genome/db_star_laevisGenome_wGTF --readFilesCommand zcat --readFilesIn /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/laevis_liver_R1_scythe.fastq.gz /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/laevis_liver_R2_scythe.fastq.gz --outFileNamePrefix /home/xue/laevis_transcriptome_mar2018/laevis_RNAseq_Star --outSAMtype BAM SortedByCoordinate
```
### Trinity: genome guided *de novo* transcriptom assembly

```
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq  --left //home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/laevis_liver_R1_scythe.fastq.gz --right /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/laevis_liver_R2_scythe.fastq.gz --genome_guided_bam /home/xue/laevis_transcriptome_mar2018/laevis_RNAseq_star/laevis_RNAseq_StarAligned.sortedByCoord.out.bam --genome_guided_max_intron 20000 --CPU 25 --full_cleanup --max_memory 200G --output /home/xue/laevis_transcriptome_mar2018/laevis_transcriptome_trinityout; echo "trinity is done at info115 in screen laevis" | mail sarahsongxy@gmail.com
```
