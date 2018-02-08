# Goal
We are building *X.laevis* transcriptome to assess the expression directionality of *X.borealis* transcripts. It means that we want to know the expression of a X.borealis transcript in X.laevis.

# Souce of RNA-seq raw data
We dig into database (ex, NCBI SRA) and literatures, the only available RNA-seq data are from Session 2016; all other are mostly RNA-seq from embryos, which doesn't have sex information. The other *X.laevis* RNA-seq are from previous sequencing of BJE4168cDNA (path to raw read: /home/xue/Xlaevis_transcriptome_QA/Laevis_TranscriptomeData/Sample_BenEvansBJE4168cDNA_Library). Those two set of raw data will be pre-processing separately before assembly and will be combined for transcriptome assembly.

# Pre-processing of Session RNA-seq
BenF downloaded raw RNA-seq data of two female liver, two female oviduct, and two male testis. Those raw data were splitted into  were trimmed with Trimmomatic, and clean with Scythe to remove 3' adaptor contamination. I have troubled trying to run Trinity with those post-processed data set. Trinity keeps giving me error messages saying that SRA data need to be clean with `SRA_TOOLKIT/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files file.sra`. I tried running the above command with the post-processed files but still got the same error message. I tried used the above commend to process the original SRA files from the database and Trinity was able to run with those. It seems like it only work with .sra files and can't handle any other file formats. Hence, I downloaded everything again, converted SRA file into fastq using sratoolkit `fast-dump` function with Trinity recommended parameter `--defline-seq '@$sn[_$rn]/$ri'`,  and clean them with Trimmomatic and Scythe:
```
#download RNAseq from SRA database; split it into R1 and R2; and processed it in proper file format that Trinity can use (time used: about 25min)
/home/xue/software/ncbi/sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump  --defline-seq '@$sn[_$rn]/$ri' --split-files SRR2515140
#gzip the downloaded file
for i in *SRR2515140*; do gzip $i; done
#trim it with trimmomatic, which was done for all through a perl script
perl run_trimmomatic.pl
#trim it again with Scythe to remove 3' contamination 
/home/xue/software/scythe-master/scythe -a /home/xue/software/scythe-master/illumina_adapters.fa -p 0.1 /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Trimmed/SRR2515154_1_R1_paired.fastq.gz | gzip > /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Scythed/SRR2515154_1_R1_paired_scythe.fastq.gz 

/home/xue/software/scythe-master/scythe -a /home/xue/software/scythe-master/illumina_adapters.fa -p 0.1 /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Trimmed/SRR2515154_2_R2_paired.fastq.gz | gzip > /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Scythed/SRR2515154_2_R2_paired_scythe.fastq.gz 

/home/xue/software/scythe-master/scythe -a /home/xue/software/scythe-master/illumina_adapters.fa -p 0.1 /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Trimmed/SRR2515140_1_R1_paired.fastq.gz | gzip > /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Scythed/SRR2515140_1_R1_paired_scythe.fastq.gz 

/home/xue/software/scythe-master/scythe -a /home/xue/software/scythe-master/illumina_adapters.fa -p 0.1 /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Trimmed/SRR2515140_2_R2_paired.fastq.gz | gzip > /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Scythed/SRR2515140_2_R2_paired_scythe.fastq.gz

#quality check with fastqc

```



# Pre-processing of BJE4168 RNA-seq
Running Trimmomatic with the same setting as what Benf did before (see Documention/Raw-reads-preparation.md). The original RNA_seq read files were seperated into small files. So I need to run trimmomatic for each subsets and I used a perl script (run_trimmomatic) to repeat the following command:  . 
```
time java -jar /home/xue/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 /home/xue/Xlaevis_transcriptome_QA/Laevis_TranscriptomeData/Sample_BenEvansBJE4168cDNA_Library/BenEvansBJE4168cDNA_Library_ACAGTG_L004_R1_001.fastq.gz /home/xue/Xlaevis_transcriptome_QA/Laevis_TranscriptomeData/Sample_BenEvansBJE4168cDNA_Library/BenEvansBJE4168cDNA_Library_ACAGTG_L004_R2_001.fastq.gz /home/xue/Xlaevis_transcriptome_QA/Laevis_TranscriptomeData/BJE4168_laevis_trimmed_data/BJE4168_L004_R1_paired.fastq.gz /home/xue/Xlaevis_transcriptome_QA/Laevis_TranscriptomeData/BJE4168_laevis_trimmed_data/BJE4168_L004_R1_unpaired.fastq.gz /home/xue/Xlaevis_transcriptome_QA/Laevis_TranscriptomeData/BJE4168_laevis_trimmed_data/BJE4168_L004_R2_paired.fastq.gz /home/xue/Xlaevis_transcriptome_QA/Laevis_TranscriptomeData/BJE4168_laevis_trimmed_data/BJE4168_L004_R2_unpaired.fastq.gz ILLUMINACLIP:/home/xue/software/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

# Building Trascriptome

