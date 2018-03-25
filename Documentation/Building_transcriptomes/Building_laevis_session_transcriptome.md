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

#quality check with fastqc
fastqc *paired*

#trim it again with Scythe to remove 3' contamination 
/home/xue/software/scythe-master/scythe -a /home/xue/software/scythe-master/illumina_adapters.fa -p 0.1 /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/laevis_session_trimmed/SRR2515140__R1_paired.fastq.gz | gzip > /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/laevis_session_scythed/SRR2515140_R1_scythed.fastq.gz 

/home/xue/software/scythe-master/scythe -a /home/xue/software/scythe-master/illumina_adapters.fa -p 0.1 /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/laevis_session_trimmed/SRR2515140__R2_paired.fastq.gz | gzip > /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/laevis_session_scythed/SRR2515140_R2_scythed.fastq.gz 


#quality check with fastqc
fastqc *
```

### Building Trascriptome
Building *De Novo* transcriptomes using Trinity with RAN-seq reads from 2 female liver samples and 2 male testis samples:
```
#with the trimmed data
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq --left /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Trimmed/*R1_paired* --right /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_Trimmed/*R2_paired* --CPU 20 --full_cleanup --max_memory 200G --output /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/Laevis_Session_TrinityOut; echo "trinity is done at info114" | mail sarahsongxy@gmail.com
#with RNA-seq for 2 female liver, 2 female oviduct, and 2 male testis
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq --left /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/laevis_session_scythed/laevis_session_R1_scythed.fastq.gz --right /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/laevis_session_scythed/laevis_session_R2_scythed.fastq.gz --CPU 25 --full_cleanup --max_memory 200G --output /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/laevis_session_trinityout; echo "trinity is done at info115" | mail sarahsongxy@gmail.com
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

### Differential expression with kallisto and EdgeR
```
#kallisto with laevis transcriptome

samplefile.txt
female    f1_liver       /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/laevis_session_scythed/SRR2515154_R1_scythed.fastq.gz /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/laevis_session_scythed/SRR2515154_R2_scythed.fastq.gz
female    f2_liver       /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/laevis_session_scythed/SRR2515140_R1_scythed.fastq.gz /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/laevis_session_scythed/SRR2515140_R2_scythed.fastq.gz
female    f1_ovary       /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/laevis_session_scythed/SRR2515157_R1_scythed.fastq.gz /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/laevis_session_scythed/SRR2515157_R2_scythed.fastq.gz
female    f2_ovary       /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/laevis_session_scythed/SRR2515143_R1_scythed.fastq.gz /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/laevis_session_scythed/SRR2515143_R2_scythed.fastq.gz
male    m1_testis      /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/laevis_session_scythed/SRR2515162_R1_scythed.fastq.gz /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/laevis_session_scythed/SRR2515162_R2_scythed.fastq.gz
male    m2_testis       /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/laevis_session_scythed/SRR2515148_R1_scythed.fastq.gz /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/laevis_session_scythed/SRR2515148_R2_scythed.fastq.gz

time perl /home/xue/software/trinityrnaseq-Trinity-v2.4.0/util/align_and_estimate_abundance.pl --transcripts /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/laevis_session_transcriptome/laevis_session_trinityout.Trinity.fasta --seqType fa --samples_file sample_file.txt --est_method kallisto --output_dir /home/xue/borealis_DE/session_laevis_deNovo_transcriptome/expression_estimate_all_kallisto --trinity_mode --prep_reference

time perl /home/xue/software/trinityrnaseq-Trinity-v2.4.0/util/abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix session2 --name_sample_by_basedir f1_liver/abundance.tsv f2_liver/abundance.tsv f1_ovary/abundance.tsv f2_ovary/abundance.tsv m1_testis/abundance.tsv m2_testis/abundance.tsv

time perl /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix session2.counts.matrix --method edgeR --samples_file sample_file.txt

time perl /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix session2.counts.matrix --method edgeR --samples_file sample_file_EdgeR.txt

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