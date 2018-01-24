# Goal
We are building *X.laevis* transcriptome to assess the expression directionality of *X.borealis* transcripts. It means that we want to know the expression of a X.borealis transcript in X.laevis.

# Souce of RNA-seq
We dig into database (ex, NCBI SRA) and literatures, the only available RNA-seq data are from Session 2016; all other are mostly RNA-seq from embryos, which doesn't have sex information. The other *X.laevis* RNA-seq are from previous sequencing of BJE4168cDNA (path to raw read: /home/xue/Xlaevis_transcriptome_QA/Laevis_TranscriptomeData/Sample_BenEvansBJE4168cDNA_Library). Those two set of raw data will be pre-processing separately before assembly and will be combined for transcriptome assembly.

# Pre-processing of Session RNA-seq

# Pre-processing of BJE4168 RNA-seq
Running Trimmomatic with the same setting as what Benf did before (see Documention/pre-processing)
```
time java -jar /home/xue/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 /home/xue/Xlaevis_transcriptome_QA/Laevis_TranscriptomeData/Sample_BenEvansBJE4168cDNA_Library/BenEvansBJE4168cDNA_Library_ACAGTG_L004_R1_001.fastq.gz /home/xue/Xlaevis_transcriptome_QA/Laevis_TranscriptomeData/Sample_BenEvansBJE4168cDNA_Library/BenEvansBJE4168cDNA_Library_ACAGTG_L004_R2_001.fastq.gz /home/xue/Xlaevis_transcriptome_QA/Laevis_TranscriptomeData/BJE4168_laevis_trimmed_data/BJE4168_L004_R1_paired.fastq.gz /home/xue/Xlaevis_transcriptome_QA/Laevis_TranscriptomeData/BJE4168_laevis_trimmed_data/BJE4168_L004_R1_unpaired.fastq.gz /home/xue/Xlaevis_transcriptome_QA/Laevis_TranscriptomeData/BJE4168_laevis_trimmed_data/BJE4168_L004_R2_paired.fastq.gz /home/xue/Xlaevis_transcriptome_QA/Laevis_TranscriptomeData/BJE4168_laevis_trimmed_data/BJE4168_L004_R2_unpaired.fastq.gz ILLUMINACLIP:/home/xue/software/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
