# RNAseq analysis
This document will use the RNAseq data from **Xenopus borealis** liver tissue of 4 males and 4 females. The samples are from a family consist of father, mother, 3 brothers, and 3 sisters. 

## Raw data
THe raw data are in fastq files.
```bash

```


## Quality Check
The quality check were done using fastqc. 
```bash
fastqc *.fastq.gz 
```
Then download files ended with `_fastqc.html` to you local computer to look at the fastqc result
```bash
scp xue@info.mcmaster.ca:/home/xue/borealis_adult_transcriptome/borealis_denovo_transcriptome_dec2018/data/trimmed/fastqc/*fastqc.html .
```


## Trimming
A perl script that was written by BenF was modified to run Trimmomatic on each inidividual rnaseq read files
```{bash}
perl run_trimmomatic.pl 
```
The Trimmomatic command in the script is:
```{bash}
java -jar /home/xue/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 $r1_name $r2_name $r1_par_out $r1_unpar_out $r2_par_out $r2_unpar_out ILLUMINACLIP:/home/xue/software/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MAXINFO:30:0.7 MINLEN:36
```
**Quality check again** We do quality check on trimmed reads again using fastqc. 
```bash
fastqc *.paired.fastq.gz 
```

## Transcriptome assembly
I merged R1 from all the samples into one files and did the same for R2 reads. 
```bash
cat *R1_paired.fastq.gz > borealis_R1_paired.fastq.gz; cat *R2_paired.fastq.gz > borealis_R2_paired.fastq.gz
```

I run trinity v2.4.0 with an additional parameter `--min_kmer_cov 2`, which should reduce the amount of low coverage transcript.
```{bash}
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq  --left /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_oct2018/Trimmed/borealis_R1_paired.fastq.gz --right /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_oct2018/Trimmed/borealis_R2_paired.fastq.gz --CPU 20 --full_cleanup --max_memory 200G --min_kmer_cov 2 --output /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_oct2018/; echo "trinity is done at info114 in screen assembly" | mail songxy2@mcmaster.ca

```
The amount of transcripts that we end up having after running the above assembly is 1.5m. The initial run in Sept 2017 had 1.6m. This Oct2018 run have 0.1m less transcripts but still very large transciptome. 

Graham (ComputeCanada) has the newest version of Trinity, hence, I build the *Xenopus borealis* transcriptome with the following bash script *build_transcriptome.sh*. 

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=1-22:00
#SBATCH --mem=210G
#SBATCH --job-name=borealis_transcriptome_dec2018
#SBATCH --account=def-ben

#build_transcriptome.sh

module load nixpkgs/16.09  gcc/7.3.0 nixpkgs/16.09
module load openmpi/3.1.2
#module load trinity/2.8.4
module load samtools/1.9
module load salmon/0.11.3
module load bowtie2/2.3.4.3
module load jellyfish/2.2.6
module load trinity/2.8.4


Trinity --seqType fq --left /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/borealis_R1_paired.fastq.gz --right /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/borealis_R2_paired.fastq.gz --CPU 20 --full_cleanup --max_memory 200G --min_kmer_cov 2 --include_supertranscripts --output /home/songxy/scratch/borealis_transcriptome_trinityOut
```
## Read mapping 
Due to lack of *X. borealis*, we mapped *X. borealis* RNAseq reads to *X. laevis* genome using STAR. 
```bash
# done by BenF
# --outFilterMismatchNmax 20 # default is 10
# --outFilterMismatchNoverLmax 0.5 # default is 0.3, ratio mismatches to mapped length
# --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 # length of reads, default is 0.66

~/bin/STAR-2.5.3a/bin/Linux_x86_64/STAR --runThreadN 10 --genomeDir /home/benf/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/STAR-Index/ --readFilesCommand zcat --outFileNamePrefix dad-star-lessStringent --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 20 --outFilterMismatchNoverLmax 0.5 --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 --readFilesIn ../Trimmed/BJE3896_dad_liver_R1_scythe.fastq.gz ../Trimmed/BJE3896_dad_liver_R2_scythe.fastq.gz
```

## Mapping transcriptome to genome

```bash
time gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92 -d laevis92_gmap -A -B 5 -t 25 -f samse --cross-species /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_dec2018/transcriptome/borealis_transcriptome_trinityOut.fasta | samtools view -S -b > /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_dec2018/analysis/transcriptome/mapping_xb_denovoTrans_xl_genomev92_gmap/borealis_denovoT_laevisv92_genome_gmap.bam
```
filter the mapping read remove the unmapped reads, which is indicated by flag 0x04. Then I piled the filtered output to bedtools, which extracts alignment coordinates for transcripts based on their CIGAR strings;
```bash
samtools view -F 0x04 -b borealis_denovoT_laevisv92_genome_gmap.bam | bedtools bamtobed -i > borealis_denovoT_laevisV92_genome_gmap_bedfile.bed
```

## Transccript expression quantification
It was done using kallisto, which is a harsh based quantification tools. 
```bash
#Indexing the transcriptomes:
kallisto index -i /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_dec2018/transcriptome/kallisto_index/borealis_transcriptome_trinityOut.fasta.kallisto_idx /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_dec2018/transcriptome/borealis_transcriptome_trinityOut.fasta

# Transcript abundance quantification:
kallisto quant -i /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta.kallisto_idx  -o female_rep7 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4082_girl_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4082_girl_liver_R2_scythe.fastq.gz)

# compile the read count for each samples into a matrix; used a script included in the Trinity RNAseq analysis package
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/util/abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix borealis_liver  --name_sample_by_basedir female_rep1/abundance.tsv female_rep2/abundance.tsv female_rep3/abundance.tsv female_rep4/abundance.tsv male_rep1/abundance.tsv male_rep2/abundance.tsv male_rep3/abundance.tsv male_rep4/abundance.tsv
```

## Differential expression analysis
There are a lot of R package that can do differential expression. I used two of them EdgeR and DESeq2. Please go to the general_scripts sub-repo for the detailed scripts.

The general steps in these two methods are:
- read in the raw transcript read counts for each sample
- filter out transcripts that have zero or very low read counts across all the samples
- put the count data into the right type of matrix 
- normalized the read counts
- fit the data into a negative binomial distribution
- perform comparative test (ex, exact test)
- extract result and store it in a output file

```bash
# DESeq2 script
bor_de_deseq2.R

# EdgeR script
bor_de_edgeR.R
```

## Post differential expression
After you have the table that contain the result of differential expression analysis, you can do the following
- filter and visualize the result
- extract the sequence of the differentially expression trancripts
```bash
#filter and visualize the result
volcano_plot.R

# extract the sequence 
perl extract_sequence.pl seq_name.file fasta.file position_of_sequence_in_seq_name.file > output

```

