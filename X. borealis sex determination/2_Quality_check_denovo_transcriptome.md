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

#### Assessing the Read Content of the Transcriptome Assembly
I ran the below script provided by Trinity to check the read contents of the transcriptome assembly. I was following the tutorial provided by Trinity (https://github.com/trinityrnaseq/trinityrnaseq/wiki/RNA-Seq-Read-Representation-by-Trinity-Assembly). Per Bowtie 2 manual, when Bowtie 2 finishes running, it prints messages summarizing what happened. These messages are printed to the “standard error” (“stderr”) filehandle, which in the below code, it was saved to file `align_stats.txt`.
```bash
#build a bowtie2 index for the transcriptome
bowtie2-build tropicalis_transcriptome_trinityOut.Trinity.fasta tropicalis_transcriptome_trinityOut.Trinity.fasta

#perform the alignment for paired-end reads 
bowtie2 -p 10 --no-unal -k 20 --phred33 -x tropicalis_transcriptome_trinityOut.Trinity.fasta -q -1 /home/xue/tropicalis_gonad_transcriptome_Dec2018/data/trim/XT_R1.fastq.gz -2 /home/xue/tropicalis_gonad_transcriptome_Dec2018/data/trim/XT_R2.fastq.gz 2>/home/xue/tropicalis_gonad_transcriptome_Dec2018/analysis/denovo_transcriptome_quality_check/read_contents_check/align_stats.txt| samtools view -@10 -Sb -o /home/xue/tropicalis_gonad_transcriptome_Dec2018/analysis/denovo_transcriptome_quality_check/read_contents_check/bowtie2.bam 
     
#capture the read alignment statistics
cat 2>&1 align_stats.txt
```
The flags used in the alignment included:
- `p/--threads NTHREADS`:Launch NTHREADS parallel search threads (default: 1).	
- `--no-unal`: Suppress SAM records for reads that failed to align.
- `-k <int>`: By default, bowtie2 searches for distinct, valid alignments for each read. When it finds a valid alignment, it continues looking for alignments that are nearly as good or better. The best alignment found is reported (randomly selected from among best if tied). Information about the best alignments is used to estimate mapping quality and to set SAM optional fields, such as AS:i and XS:i. When -k is specified, however, bowtie2 behaves differently. Instead, it searches for at most <int> distinct, valid alignments for each read. The search terminates when it can’t find more distinct valid alignments, or when it finds <int>, whichever happens first.
- `--phred33`: Input qualities are ASCII chars equal to the Phred quality plus 33. This is also called the “Phred+33” encoding, which is used by the very latest Illumina pipelines.
- `-x <bt2-idx>`: The basename of the index for the reference genome.
- `-q`: Reads are FASTQ files.
- `-1 <m1>`: Comma-separated list of files containing mate 1s
- `-2 <m2>`: Comma-separated list of files containing mate 2s


Below is the result:
```

```
Note/interpretation about the result
- good
- need to be cautious about
