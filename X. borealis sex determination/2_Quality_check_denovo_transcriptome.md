# Quality Assessment for *de novo* Assembled Transcriptome 
## Transfer Data to Info
```bash
#I transfer the transcriptome and supertranscriptome data, and mapping result back to Info for down stream analysis
#since on Info, jobs will be ran more immediately without the queueing system, which Graham use and sometime jobs just queueing for a day or more.
scp -r songxy@graham.computecanada.ca:/home/songxy/projects/def-ben/songxy/tropicalis_gonad_transcriptome/data/tropicalis_gonad_supertranscriptome_dec2018/ /home/xue/tropicalis_gonad_transcriptome_Dec2018/data/tropicali_gonad_transcriptome_trinityOut

scp -r songxy@graham.computecanada.ca:/home/songxy/projects/def-ben/songxy/tropicalis_gonad_transcriptome/data/tropicalis_transcriptome_build_dec2018 /home/xue/tropicalis_gonad_transcriptome_Dec2018/data/tropicali_gonad_transcriptome_trinityOut
```
## Basic check 
I followed the Trinity instruction and ran a few protocols provided by Trinity. 

## Assessing the Read Content of the Transcriptome Assembly
I ran the below script provided by Trinity to check the read contents of the transcriptome assembly. I was following the tutorial provided by Trinity (https://github.com/trinityrnaseq/trinityrnaseq/wiki/RNA-Seq-Read-Representation-by-Trinity-Assembly).
```bash
#build a bowtie2 index for the transcriptome
bowtie2-build tropicalis_transcriptome_trinityOut.Trinity.fasta tropicalis_transcriptome_trinityOut.Trinity.fasta

#perform the alignment (example for paired-end reads) 
bowtie2 -p 10 -q --no-unal -k 20 -x tropicalis_transcriptome_trinityOut.Trinity.fasta -1 /home/xue/tropicalis_gonad_transcriptome_Dec2018/data/trim/XT_R1.fastq.gz -2 /home/xue/tropicalis_gonad_transcriptome_Dec2018/data/trim/XT_R2.fastq.gz 2 > align_stats.txt| samtools view -@10 -Sb -o bowtie2.bam 
     
#capture the read alignment statistics
cat 2>&1 align_stats.txt
```
Below is the result:
```

```
Note/interpretation about the result
     - good
     - need to be cautious about
