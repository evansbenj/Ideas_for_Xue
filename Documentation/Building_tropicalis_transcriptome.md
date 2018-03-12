# Goal
We are building *X.laevis* transcriptome to assess the expression directionality of *X.borealis* transcripts. It means that we want to know the expression of a X.borealis transcript in X.laevis.

# General pipeline
The general pipeline is very similar for processing *X. tropicalis* and *X. laevis* data.
- Obtaining raw data
- Trimming with Trimmomatic and Scyth; assessing with fastqc
- Transcriptome assemble with Trinity
- Expression level quantification with Kallisto
- Convert to TPM with R-script
- Differential expression analysis with EdgeR

# *X. laevis* with Session data
### Souce of RNA-seq raw data
We dig into database (ex, NCBI SRA) and literatures, the only available RNA-seq data are from Session 2016; all other are mostly RNA-seq from embryos, which doesn't have sex information. The other *X.laevis* RNA-seq are from previous sequencing of BJE4168cDNA (path to raw read: /home/xue/Xlaevis_transcriptome_QA/Laevis_TranscriptomeData/Sample_BenEvansBJE4168cDNA_Library). Those two set of raw data will be pre-processing separately before assembly and will be combined for transcriptome assembly.

### Pre-processing of Session RNA-seq
BenF downloaded raw RNA-seq data of two female liver, two female oviduct, and two male testis. Those raw data were splitted into  were trimmed with Trimmomatic, and clean with Scythe to remove 3' adaptor contamination. I have troubled trying to run Trinity with those post-processed data set. Trinity keeps giving me error messages saying that SRA data need to be clean with `SRA_TOOLKIT/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files file.sra`. I tried running the above command with the post-processed files but still got the same error message. I tried used the above commend to process the original SRA files from the database and Trinity was able to run with those. It seems like it only work with .sra files and can't handle any other file formats. Hence, I downloaded everything again, converted SRA file into fastq using sratoolkit `fast-dump` function with Trinity recommended parameter `--defline-seq '@$sn[_$rn]/$ri'`,  and clean them with Trimmomatic and Scythe:
```
#download RNAseq from SRA database; split it into R1 and R2; and processed it in proper file format that Trinity can use (time used: about 25min)
/home/xue/software/ncbi/sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump  --defline-seq '@$sn[_$rn]/$ri' --split-files --gzip SRR2515140
/home/xue/software/ncbi/sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump  --defline-seq '@$sn[_$rn]/$ri' --split-files --gzip SRR2515154

#trim them with trimmomatic, which was done for all through a perl script
perl run_trimmomatic.pl
#trim each of them individually  with trimmomatic (about 35min)
time java -jar /home/xue/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ../Laevis_Session_RawData/SRR2515154_1.fastq.gz ../Laevis_Session_RawData/SRR2515154_1.fastq.gz SRR2515154__R1_paired.fastq.gz SRR2515154__R1_unpaired.fastq.gz SRR2515154__R2_paired.fastq.gz SRR2515154__R2_unpaired.fastq.gz ILLUMINACLIP:/home/xue/software/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#trim it again with Scythe to remove 3' contamination 
/home/xue/software/scythe-master/scythe -a /home/xue/software/scythe-master/illumina_adapters.fa -p 0.1 /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Trimmed/SRR2515154_1_R1_paired.fastq.gz | gzip > /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Scythed/SRR2515154_1_R1_paired_scythe.fastq.gz 

/home/xue/software/scythe-master/scythe -a /home/xue/software/scythe-master/illumina_adapters.fa -p 0.1 /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Trimmed/SRR2515154_2_R2_paired.fastq.gz | gzip > /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Scythed/SRR2515154_2_R2_paired_scythe.fastq.gz 

/home/xue/software/scythe-master/scythe -a /home/xue/software/scythe-master/illumina_adapters.fa -p 0.1 /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Trimmed/SRR2515140_1_R1_paired.fastq.gz | gzip > /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Scythed/SRR2515140_1_R1_paired_scythe.fastq.gz 

/home/xue/software/scythe-master/scythe -a /home/xue/software/scythe-master/illumina_adapters.fa -p 0.1 /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Trimmed/SRR2515140_2_R2_paired.fastq.gz | gzip > /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Scythed/SRR2515140_2_R2_paired_scythe.fastq.gz

/home/xue/software/scythe-master/scythe -a /home/xue/software/scythe-master/illumina_adapters.fa -p 0.1 /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Trimmed/SRR2515154_2_R2_unpaired.fastq.gz | gzip > /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Scythed/SRR2515154_2_R2_unpaired_scythe.fastq.gz 

#quality check with fastqc
fastqc *paired*
```

### Building Trascriptome
Building *De Novo* transcriptomes using Trinity with RAN-seq reads from 2 female liver samples and 2 male testis samples:
```
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq --left /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Trimmed/*R1_paired* --right /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Trimmed/*R2_paired* --CPU 20 --full_cleanup --max_memory 200G --output /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_TrinityOut; echo "trinity is done at info114" | mail sarahsongxy@gmail.com

```
Build individual transcriptome for each sample
```
#liver female 1
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq --left /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Trimmed/SRR2515140__R1_paired.fastq.gz --right /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Trimmed/SRR2515140__R2_paired.fastq.gz --CPU 20 --full_cleanup --max_memory 200G --output /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/laevis_session_transcriptome/individual/liver_f1_trinityout; echo "trinity is done in screen liver at info114" | mail sarahsongxy@gmail.com
#liver female 2
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq --left /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Trimmed/SRR2515154__R1_paired.fastq.gz --right /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Trimmed/SRR2515154__R2_paired.fastq.gz --CPU 20 --full_cleanup --max_memory 200G --output /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/laevis_session_transcriptome/individual/liver_f2_trinityout; echo "trinity is done in screen liver at info115" | mail sarahsongxy@gmail.com
#testis male 1
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq --left /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Trimmed/SRR2515162__R1_paired.fastq.gz --right /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Trimmed/SRR2515162__R2_paired.fastq.gz --CPU 20 --full_cleanup --max_memory 200G --output /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/laevis_session_transcriptome/individual/testis_m1_trinityout; echo "trinity is done in screen testis at info114" | mail sarahsongxy@gmail.com
#testis male 2
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq --left /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Trimmed/SRR2515148__R1_paired.fastq.gz --right /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Trimmed/SRR2515148__R2_paired.fastq.gz --CPU 20 --full_cleanup --max_memory 200G --output /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/laevis_session_transcriptome/individual/testis_m2_trinityout; echo "trinity is done in screen testis at info115" | mail sarahsongxy@gmail.com
```

### Pre-processing of BJE4168 RNA-seq (Abadoned, unknow sex)
Running Trimmomatic with the same setting as what Benf did before (see Documention/Raw-reads-preparation.md). The original RNA_seq read files were seperated into small files. So I need to run trimmomatic for each subsets and I used a perl script (run_trimmomatic.pl) to repeat the following command:  . 
```
time java -jar /home/xue/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 /home/xue/Xlaevis_transcriptome_QA/Laevis_TranscriptomeData/Sample_BenEvansBJE4168cDNA_Library/BenEvansBJE4168cDNA_Library_ACAGTG_L004_R1_001.fastq.gz /home/xue/Xlaevis_transcriptome_QA/Laevis_TranscriptomeData/Sample_BenEvansBJE4168cDNA_Library/BenEvansBJE4168cDNA_Library_ACAGTG_L004_R2_001.fastq.gz /home/xue/Xlaevis_transcriptome_QA/Laevis_TranscriptomeData/BJE4168_laevis_trimmed_data/BJE4168_L004_R1_paired.fastq.gz /home/xue/Xlaevis_transcriptome_QA/Laevis_TranscriptomeData/BJE4168_laevis_trimmed_data/BJE4168_L004_R1_unpaired.fastq.gz /home/xue/Xlaevis_transcriptome_QA/Laevis_TranscriptomeData/BJE4168_laevis_trimmed_data/BJE4168_L004_R2_paired.fastq.gz /home/xue/Xlaevis_transcriptome_QA/Laevis_TranscriptomeData/BJE4168_laevis_trimmed_data/BJE4168_L004_R2_unpaired.fastq.gz ILLUMINACLIP:/home/xue/software/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
Running Trinity:
```
/home/xue/software/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq --left /home/xue/Trinity_Run2/BJE4168_laevis_trimmed_data/BJE4168_R1_paired_fastq.gz --right /home/xue/Trinity_Run2/BJE4168_laevis_trimmed_data/BJE4168_R2_paired_fastq.gz --CPU 20 --full_cleanup --max_memory 200G --output /home/xue/Trinity_Run2/BJE4168_laevis_trimmed_data
```




# Tropicalis with Marin data
We are using the RNA-seq data of X. tropicalis as a outgroup for an analysis where we examine directionality of expression level (up-regulated or down-regulated) in X. borealis. 

### RNA-seq raw data from SRA
here are RNA-seq data that are from adult liver tissues and have sex information but the researchers did single end sequencing (BioProject:PRJNA381064). Those are from a study about dosage compensation mechanism in a reptile lineage. (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5741051/). We download RNA-seq data for 2 female liver samples (SRR5412274 and SRR5412273) and 2 male liver samples (SRR5412276 and SRR5412275). 
```
# run bash script to download 4 samples
. download_sra_rnaseq.sh
#command(in the script) to download SRA files  

/home/xue/software/ncbi/sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump  --defline-seq '@$sn[_$rn]/$ri' --gzip SRR5412274
```
### Transcriptome
Building transcriptomes with 2 liver RNA-seq data:
```
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq --samples_file /home/xue/borealis_DE/tropicalis_transcriptome/tropicalis_trimmed/tropicalis_samples_file --CPU 20 --full_cleanup --max_memory 200G --output /home/xue/borealis_DE/tropicalis_transcriptome/tropicalis_trinityout; echo "trinity is done at info114" | mail sarahsongxy@gmail.com

```





