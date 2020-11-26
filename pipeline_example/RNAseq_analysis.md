# RNAseq analysis
This document will use the RNAseq data from *Xenopus borealis* liver tissue of 4 males and 4 females. The samples are from a family consist of father, mother, 3 brothers, and 3 sisters. 

## Raw data
THe raw data are in fastq files.
```bash
#path to the raw reads
/home/benf/Borealis-Family-Transcriptomes-July2017/Data/RawReads
```
One sequence read from dad's liver R1 file

- line 1: name of read; 
- line 2: sequence of the read
- line 3: a plus sign
- line 4: quality score of the base calling from the sequencer; this tells you how confident the sequencer is calling when it call the base A/T/C/G. 
```
@D00430:256:CB2H9ANXX:1:1106:1227:2076 1:N:0:ATCACG
GGACTCAGACTGATTCCTGCACAGGGGGTTCATACCCATCAGCACGTTCCACTTTGGCACCAGTCTCGTCTCCTGCTGGTTTTGCAGTACCACCACCCTCGCCATGAAGTTCCATTAGTTTGCCCA
+
BCBBCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG0<FGEGGGGGGGGGGGGGGGGGGGGGGGDGGGGGGEGEGFGGGEGGGGGGBGGGGGGGGGGGGFECGGGGCEGGGEGGG/
```

## Quality Check
The quality check were done using fastqc. 
```bash
fastqc *.fastq.gz 
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
Then download files ended with `_fastqc.html` to my local computer to look at the fastqc result.

-usually RNAseq data have high duplication level and it's not a big issue to worry about even when the duplication session gives a error message 

```bash
#download fastqc result to local computer
scp xue@info.mcmaster.ca:/home/xue/borealis_adult_transcriptome/borealis_denovo_transcriptome_dec2018/data/trimmed/fastqc/*fastqc.html .
```

## Read mapping 
Due to lack of *X. borealis* genome, we mapped *X. borealis* RNAseq reads to *X. laevis* genome using STAR. 
```bash
# done by BenF
# --outFilterMismatchNmax 20 # default is 10
# --outFilterMismatchNoverLmax 0.5 # default is 0.3, ratio mismatches to mapped length
# --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 # length of reads, default is 0.66

~/bin/STAR-2.5.3a/bin/Linux_x86_64/STAR --runThreadN 10 --genomeDir /home/benf/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/STAR-Index/ --readFilesCommand zcat --outFileNamePrefix dad-star-lessStringent --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 20 --outFilterMismatchNoverLmax 0.5 --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 --readFilesIn ../Trimmed/BJE3896_dad_liver_R1_scythe.fastq.gz ../Trimmed/BJE3896_dad_liver_R2_scythe.fastq.gz
```

## Transcriptome assembly
I merged R1 reads from all inividual into one files and did the same for R2 reads. 
```bash
cat *R1_paired.fastq.gz > borealis_R1_paired.fastq.gz; cat *R2_paired.fastq.gz > borealis_R2_paired.fastq.gz
```

I run trinity v2.4.0 with an additional parameter `--min_kmer_cov 2`, which should reduce the amount of low coverage transcript.
```{bash}
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq  --left /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_oct2018/Trimmed/borealis_R1_paired.fastq.gz --right /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_oct2018/Trimmed/borealis_R2_paired.fastq.gz --CPU 20 --full_cleanup --max_memory 200G --min_kmer_cov 2 --output /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_oct2018/; echo "trinity is done at info114 in screen assembly" | mail songxy2@mcmaster.ca

```
The amount of transcripts that we end up having after running the above assembly (OCt2018) is 1.5m. The initial run in Sept 2017 had 1.6m. This Oct2018 run have 0.1m less transcripts but still very large transciptome. 

Graham (ComputeCanada) has the newest version of Trinity, hence, I build the *Xenopus borealis* transcriptome with the following bash script *build_transcriptome.sh*. 

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=3-22:00
#SBATCH --mem=260G
#SBATCH --job-name=borealis_transcriptome_Jun2020
#SBATCH --account=def-ben

module load nixpkgs/16.09  gcc/7.3.0 nixpkgs/16.09
module load openmpi/3.1.2
module load samtools/1.9
module load salmon/0.11.3
module load bowtie2/2.3.4.3
module load jellyfish/2.2.6
module load trinity/2.8.4


Trinity --seqType fq --left /home/songxy/scratch/borealis_transcriptome/trimmed/borealis_liver_R1_paired.fastq.gz --right /home/songxy/scratch/borealis_transcriptome/trimmed/borealis_liver_R2_paired.fastq.gz --CPU 20 --full_cleanup --max_memory 250G --min_kmer_cov 2 --output /home/songxy/scratch/borealis_transcriptome/borealis_adult_liver_transcriptome_trinityOut
```


## Mapping transcriptome to genome
Since genome contains introns+exons and transcriptome only contain exons, we need a gap-aware mapper such [Gmap](https://academic.oup.com/bioinformatics/article/21/9/1859/409207) to align *X. borealis* transcripts to *X. laevis* genome.  
```bash
# index the X. laevis genome
# default k-mer value for genomic index is 15; use -k to change it to int <= 16
# gmap_build [options...] -d <genomename> <fasta_files>
gmap_build -d db_gmap_xl92 /home/xue/genome_data/laevis_genome/XL9_2.fa 

# map X. borealis transcript to X. laevis genome
# added --cross-species to use a more sensitive search for canonical splicing
# piled to samtools to output as a bamfile to save space on cluster
# ran on info115; this took about 37hrs
# real    2269m54.858s
# user    32242m29.231s
# sys     21379m39.735s
time gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92 -d laevis92_gmap -A -B 5 -t 25 -f samse --cross-species /home/xue/borealis_adult_transcriptome/borealis_liver_denovo_transcriptome_May2020/build_transcriptome/borealis_adult_liver_transcriptome_trinityOut.Trinity.fasta | samtools view -S -b > /home/xue/borealis_adult_transcriptome/borealis_liver_denovo_transcriptome_May2020/mapping_xb_liverTrans_laevisG92_gmap/borealis_liverTrans_laevisGenomeV92_gmap.bam
```
filtering alignment: 

1) remove the unmapped transcripts, which is indicated by flag 0x04. I piled the filtered output to bedtools, which extracts alignment coordinates for transcripts based on their CIGAR strings.
2) filter the bedfile and keep only the hit with the highest mapping quality score; if there are two hit with the same mapping quality score and mapped to the same chromosome, keep the first one only;

```{bash}
#1) on info: removed those that didn't mapped
samtools view -F 0x04 -b borealis_liverTrans_borGenomeV1_gmap.bam | bedtools bamtobed -i > borealis_liverTrans_borGenomeV1_gmap.bed

#2) run the below script to keep the top hit
Rscript ~/script/filter_bedfile_keep1match.R borealis_liverTrans_laevisGenomeV92_gmap.bed borealis_liverTrans_laevisGenomeV92_gmap_keep1match.tsv
```

## Transccript expression quantification
It was done using kallisto. 

- Kallisto is a RNA-seq transcript abundance quantification tool that estimates transcripts abundance with pseudo-alignment method
- Input: transcriptome file and trimmed reads
- output: plain text file containing transcript abundances information for each transcript - transcript id, length, effective length, estimated count, count in TPM (transcript in millions); the output file is used in downstream analysis such as differential expression level analysis
- below is a brief description of how it works in steps:
  - indexing step: build transcriptome de Bruijn Graph (T-DBG) with k-mers from transcripts in the assembled transcriptome. In a hash table, store the information about the k-mers - which transcripts a k-mer is from and along with information the position of the k-mer in the transcript. 
  - pseudo alignment step: for each k-mer in the read, look for those k-mers in the hash table to found out which transcripts those k-mers could be from. In the case where a read has 5 k-mers and k-mer #1-#3 are from transcript A & B & C but k-mer #4-#5 are from transcript A & B, kallisto will then decide that the read is probably from transcript A &B.      
  - quantification step:  for a multi-mapping read like the example above, the read is fractionally assigned to each transcript based on the expectation-maximization algorithm. With the above example, a fraction of the read will be assigned to transcript A and transcript B.

```bash
# Indexing the transcriptome:
kallisto index -i /home/xue/borealis_adult_transcriptome/borealis_liver_denovo_transcriptome_May2020/raw_count_kallisto/borealis_adult_liver_transcriptome_trinityout.fasta.kallisto_idx /home/xue/borealis_adult_transcriptome/borealis_liver_denovo_transcriptome_May2020/build_transcriptome/borealis_adult_liver_transcriptome_trinityout.fasta

# Transcript abundance quantification:
kallisto quant -i /home/xue/borealis_adult_transcriptome/borealis_liver_denovo_transcriptome_May2020/raw_count_kallisto/borealis_adult_liver_transcriptome_trinityout.fasta.kallisto_idx  -o BJE4039_boy <(gunzip -c /home/xue/borealis_adult_transcriptome/borealis_denovo_transcriptome_dec2018/data/trimmed/BJE4039_boy_liver_R1_paired.fastq.gz) <(gunzip -c /home/xue/borealis_adult_transcriptome/borealis_denovo_transcriptome_dec2018/data/trimmed/BJE4039_boy_liver_R2_paired.fastq.gz)

# compile the read count for each samples into a matrix; used a script included in the Trinity RNAseq analysis package
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/util/abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix borealis_liver  --name_sample_by_basedir BJE3897_mom/abundance.tsv BJE4009_girl/abundance.tsv BJE4072_girl/abundance.tsv BJE4082_girl/abundance.tsv BJE3896_dad/abundance.tsv BJE3929_boy/abundance.tsv BJE4017_boy/abundance.tsv BJE4039_boy/abundance.tsv
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
The critiria used to identify significantlly differential expressed transcripts are:

- FDR or adjusted_p_value <= 0.05
- logFC >= 2 and logFC <= -2

## Post differential expression analysis
After you have the table that contain the result of differential expression analysis, you can do the following
- filter and visualize the result
- extract the sequence of the differentially expression trancripts
```bash
#filter and visualize the result
volcano_plot.R

# extract the sequence 
perl extract_sequence.pl seq_name.file fasta.file position_of_sequence_in_seq_name.file > output
```
## Links to useful articles
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6373869/
- https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
- https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8


