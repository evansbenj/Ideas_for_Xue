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


# *X. Tropicalis* with Marin et al data
We are using the RNA-seq data of X. tropicalis as a outgroup for an analysis where we examine directionality of expression level (up-regulated or down-regulated) in X. borealis. 

### RNA-seq raw data from SRA + Pre-processing
here are RNA-seq data that are from adult liver tissues and have sex information but the researchers did single end sequencing (BioProject:PRJNA381064). Those are from a study about dosage compensation mechanism in a reptile lineage. (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5741051/). We download RNA-seq data for 2 female liver samples (SRR5412274 and SRR5412273) and 2 male liver samples (SRR5412276 and SRR5412275). 
```
# run bash script to download 4 samples
sra_single_trimmomatic.pl
#commands(in the script) to download SRA files  
/home/xue/software/ncbi/sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump  --defline-seq '@$sn[_$rn]/$ri' --gzip SRR5412274
#trimmomatic command for individual sample in the script
time java -jar /home/xue/software/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 $infile_name $outfile_name ILLUMINACLIP:/home/xue/software/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
### Building tropicalis transcriptome using Trinity (*de novo*)
Building *de novo* transcriptomes with 4 liver RNA-seq data:
```
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq --samples_file /home/xue/borealis_DE/tropicalis_transcriptome/tropicalis_trimmed/tropicalis_samples_file --CPU 20 --full_cleanup --max_memory 200G --output /home/xue/borealis_DE/tropicalis_transcriptome/tropicalis_trinityout; echo "trinity is done at info114" | mail sarahsongxy@gmail.com

```
### Building tropicalis transcriptome using Trinity (genome guided)
Map RNAseq to tropicalis genome (v9.1). 
```
#getting tropicalis genome v9.1, gff3 and GTF files
wget ftp://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/XT9_1.fa.gz
wget ftp://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/XENTR_9.1_Xenbase.GTF
wget ftp://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/XENTR_9.1_Xenbase.gff3

#building Star Index
STAR --runThreadN 20 --runMode genomeGenerate --genomeDir /home/xue/genome_data/tropicalis_genome/db_tropicalis_star --genomeFastaFiles /home/xue/genome_data/tropicalis_genome/XT9_1.fa --sjdbGTFfile  /home/xue/genome_data/tropicalis_genome/XENTR_9.1_Xenbase.GTF   --genomeChrBinNbits 16

#mapping tropicalis RNAseq to tropicalis genome to get the SortedbyCoordinate.bam file
STAR --runThreadN 20 --genomeDir ~/genome_data/tropicalis_genome/db_tropicalis_star --readFilesCommand zcat --readFilesIn /home/xue/tropicalis_transcriptome/tropicalis_scythed/tropicalis_sra_scythed.fastq.gz  --outFileNamePrefix /home/xue/tropicalis_transcriptome/tropicalis_RNAseq_genome_Star/tropicalis_mapping_RNAseq_genome --outSAMtype BAM SortedByCoordinate

```
Building *de novo* transcriptomes with 4 liver RNA-seq data:


### Differential expression analysis
Kallisto
```
#building kallisto index
kallisto index -i /home/xue/borealis_DE/tropicalis_transcriptome/tropicalis_kallisto/tropicalis_trinityout.idx /home/xue/borealis_DE/tropicalis_transcriptome/tropicalis_trinityout.Trinity.fasta

#kallisto quantification for each individual file
kallisto quant -i /home/xue/borealis_DE/tropicalis_transcriptome/tropicalis_kallisto/tropicalis_trinityout.idx -o /home/xue/borealis_DE/tropicalis_transcriptome/tropicalis_kallisto/SRR5412273 -b 100 --single -l 100 -s 20 <(gunzip -c /home/xue/borealis_DE/tropicalis_transcriptome/tropicalis_trimmed/SRR5412273_trimmed.fastq.gz)
kallisto quant -i /home/xue/borealis_DE/tropicalis_transcriptome/tropicalis_kallisto/tropicalis_trinityout.idx -o /home/xue/borealis_DE/tropicalis_transcriptome/tropicalis_kallisto/SRR5412274 -b 100 --single -l 100 -s 20 <(gunzip -c /home/xue/borealis_DE/tropicalis_transcriptome/tropicalis_trimmed/SRR5412274_trimmed.fastq.gz)
kallisto quant -i /home/xue/borealis_DE/tropicalis_transcriptome/tropicalis_kallisto/tropicalis_trinityout.idx -o /home/xue/borealis_DE/tropicalis_transcriptome/tropicalis_kallisto/SRR5412275 -b 100 --single -l 100 -s 20 <(gunzip -c /home/xue/borealis_DE/tropicalis_transcriptome/tropicalis_trimmed/SRR5412275_trimmed.fastq.gz)
kallisto quant -i /home/xue/borealis_DE/tropicalis_transcriptome/tropicalis_kallisto/tropicalis_trinityout.idx -o /home/xue/borealis_DE/tropicalis_transcriptome/tropicalis_kallisto/SRR5412276 -b 100 --single -l 100 -s 20 <(gunzip -c /home/xue/borealis_DE/tropicalis_transcriptome/tropicalis_trimmed/SRR5412276_trimmed.fastq.gz)


Or

#running kallisto using Trinity script
time perl /home/xue/software/trinityrnaseq-Trinity-v2.4.0/util/align_and_estimate_abundance.pl --transcripts /home/xue/borealis_DE/tropicalis_transcriptome/tropicalis_trinityout.Trinity.fasta --seqType fa --samples_file /home/xue/borealis_DE/tropicalis_transcriptome/tropicalis_kallisto/tropicalis_samplefile.txt --est_method kallisto --output_dir /home/xue/borealis_DE/tropicalis_transcriptome/tropicalis_kallisto --output_prefix tropicalis --trinity_mode --prep_reference
```





