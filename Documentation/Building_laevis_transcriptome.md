

# Building laevis transcriptome with Trinity (genome-guided method)

### Pre-assembly: mapping reads to laevis genome
```
STAR --runThreadN 20 --genomeDir ~/genome_data/laevis_genome/db_star_laevisGenome_wGTF --readFilesCommand zcat --readFilesIn /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/laevis_liver_R1_scythe.fastq.gz /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/laevis_liver_R2_scythe.fastq.gz --outFileNamePrefix /home/xue/laevis_transcriptome_mar2018/laevis_RNAseq_Star --outSAMtype BAM SortedByCoordinate
```
### Trinity: genome guided *de novo* transcriptom assembly
```
```
