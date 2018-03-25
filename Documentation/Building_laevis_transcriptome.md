

# Building laevis transcriptome with Trinity (genome-guided method)

###Pre-assembly: mapping reads to laevis genome
```
STAR --runThreadN 20 --genomeDir /home/xue/borealis_DE/xl_genome/STAR_XLgenome \ 
--readFilesIn /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/laevis_liver_R1_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/laevis_liver_R2_scythe.fastq.gz  \
--genomeLoad LoadAndKeep --outFileNamePrefix ./xbxl_ --outSAMtype BAM Unsorted
```
