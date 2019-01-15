
## Obtaining RNAseq raw data
### Downloading Raw RNAseq data from sequencing center
```
```
### Upload the Raw RNAseq data to info
```
```
### Do fastqc on raw data
```
```

## Scythe on raw data
more detail about scythe: https://github.com/vsbuffalo/scythe. 

#### how to run Scythe
```
/home/xue/software/scythe-master/scythe -a /home/xue/software/scythe-master/illumina_adapters.fa -p 0.1 -o /home/xue/tropicalis_gonad_transcriptome_Dec2018/raw_data/XT1_R1_scythe.fastq.gz /home/xue/tropicalis_gonad_transcriptome_Dec2018/raw_data/XT1_R1.fastq.gz
```
What are the parameters
- -a: adapter file in fasta format
- -p: the prior contamination rate is 0.05; BenF used 0.1.
- -o: output file name

Other options:
- -q: Illumina's quality scheme (pipeline > 1.3) is used. Sanger or Solexa (pipeline < 1.3) qualities can be specified
- -n: one can specify the minimum match length argument with -n <integer> and 
- -M: the minimum length of sequence (discarded less than or equal to this parameter) to keep after trimming with -M <integer>

#### how to run Scythe on all the file with one command
```
for  in *fastq.gz ; do name=$(grep -o "XT[0-9]*_[A-Z][0-9]" <(echo $i)); /home/xue/software/scythe-master/scythe -a /home/xue/software/scythe-master/illumina_adapters.fa -p 0.1 $i | gzip > /home/xue/tropicalis_gonad_transcriptome_Dec2018/scythed_data/$name\_scythe.fastq.gz; done
```
#### run fastqc on scythe data
We did it for all fastq.gz file in the folder
```
for i in *fastq.gz ; do fastqc $i; done
```

and then download the file to a local folder and then inspect with a browser like this:

```
scp xue@info.mcmaster.ca:/home/xue/tropicalis_gonad_transcriptome_Dec2018/scythed_data/*html .
```

## Trimmomatic
run it like this if it is paired end
learn how to run trimmomatic and what the command means: http://www.usadellab.org/cms/?page=trimmomatic

```bash
java -jar /home/xue/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -trimlog trim_out.txt /home/xue/transcriptome_data/Sample_BenEvansBJE3909cDNA_Library/BenEvansBJE3909cDNA_Library_GTGAAA_L004_R1_001.fastq.gz /home/xue/transcriptome_data/Sample_BenEvansBJE3909cDNA_Library/BenEvansBJE3909cDNA_Library_GTGAAA_L004_R2_001.fastq.gz /home/xue/transcriptome_data/BJE3909_tropicalis_trimmed_data/BenEvansBJE3909cDNA_Library_GTGAAA_L004_R1_001_paired.fastq.gz /home/xue/transcriptome_data/BJE3909_tropicalis_trimmed_data/BenEvansBJE3909cDNA_Library_GTGAAA_L004_R1_001_single.fastq.gz /home/xue/transcriptome_data/BJE3909_tropicalis_trimmed_data/BenEvansBJE3909cDNA_Library_GTGAAA_L004_R2_001_paired.fastq.gz /home/xue/transcriptome_data/BJE3909_tropicalis_trimmed_data/BenEvansBJE3909cDNA_Library_GTGAAA_L004_R2_001_single.fastq.gz ILLUMINACLIP:/home/xue/Trimmomatic-0.36/adapters/TruSeq2-PE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36
```
What the parameters mean:
- LLUMINACLIP: Cut adapter and other illumina-specific sequences from the read.
- SLIDINGWINDOW: Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold.
- LEADING: Cut bases off the start of a read, if below a threshold quality
- TRAILING: Cut bases off the end of a read, if below a threshold quality
- CROP: Cut the read to a specified length
- HEADCROP: Cut the specified number of bases from the start of the read
- MINLEN: Drop the read if it is below a specified length
- TOPHRED33: Convert quality scores to Phred-33
- TOPHRED64: Convert quality scores to Phred-64

#### run fastqc on scythe data
We did it for all fastq.gz file in the folder
```bash
for i in *fastq.gz ; do fastqc $i; done
```


## Transcriptome assembly with Trinity
More detail about trinity here: https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity

Before I start Trinity run, I will like to combine all the R1_paired together and all the R2_paired together by doing:
```bash
cat *_R1_paired.fastq.gz > XT_R1.fastq.gz
cat *_R2_paired.fastq.gz > XT_R2.fastq.gz
```

Since I can install Trinity 2.8 and 2.7 after multiple tries, I decide to go with version 2.5.1 on info. 
```bash
export PATH="/home/xue/software/bowtie2-2.3.4.3:$PATH";time /home/xue/software/trinityrnaseq-Trinity-v2.5.1/Trinity --seqType fq  --left /home/xue/tropicalis_gonad_transcriptome_Dec2018/trim/XT_R1.fastq.gz --right /home/xue/tropicalis_gonad_transcriptome_Dec2018/trim/XT_R2.fastq.gz --CPU 20 --inchworm_cpu 6 --full_cleanup --max_memory 200G --min_kmer_cov 2 --output /home/xue/tropicalis_gonad_transcriptome_Dec2018/tropicali_gonad_transcriptome_trinityOut;echo "trinity done in screen tropicalis_gonad on info114"|mail songxy2@mcmaster.ca
```


Graham alread have the newest version of Trinity (Trinity 2.8.4) installed, I am running the newest verion on Graham with the below bash script.
```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=5-12:00
#SBATCH --mem=210G
#SBATCH --job-name=tropicalis_transcriptome_dec2018
#SBATCH --account=def-ben

module load nixpkgs/16.09  gcc/7.3.0 nixpkgs/16.09
module load openmpi/3.1.2
#module load trinity/2.8.4
module load samtools/1.9
module load salmon/0.11.3
module load bowtie2/2.3.4.3
module load jellyfish/2.2.6
module load trinity/2.8.4

Trinity --seqType fq --left /home/songxy/projects/def-ben/songxy/tropicalis_gonad_transcriptome/data/trim/XT_R1.fastq.gz --right /home/songxy/projects/def-ben/songxy/tropicalis_gonad_transcriptome/data/trim/XT_R2.fastq.gz --CPU 20 --full_cleanup --max_memory 200G --min_kmer_cov 2 --include_supertranscripts --output /home/songxy/scratch/tropicalis_transcriptome/tropicalis_transcriptome_trinityOut
```
The one ran in info gave my error message and stopped running. The one in graham ran for a little bit more than 1 days with max 191Gb of memories. Much shorter time than I expected. 

### output from trinity
- the output will be a fasta file by default. It will have the transcript name and transcript sequence for each transcript it builds.
- example of a transcript name is: >TRINITY_DN86389_c0_g1_i1 len=398 path=[0:0-159 2:160-208 3:209-275 5:276-397]
  - the c means cluster
  - the g means gene
  - the i means isoform: different isoform of a gene are supposely splice variants if Trinity did the assembly correctly.

# Mapping to *X.tropicalis* genome (v91)
I mapped the *de novo* assembled tropicalis transcriptome to the tropicalis genome (v91) using gmap on Graham. 
```bash
#downloading the genome from Xenbase (ftp://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/) 
wget ftp://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.1/XT9_1.fa.gz

#indexing
gmap_build -d db_gmap_tropicalis_v91/ -D /home/songxy/scratch/tropicalis_transcriptome/tropicalis_genome/db_gmap_tropicalis_v91 -g /home/songxy/scratch/tropicalis_transcriptome/tropicalis_genome/XT9_1.fa.gz

#mapping 
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --time=01:20:00
#SBATCH --mem=10G
#SBATCH --job-name=tropicalis_gmap_indexing
#SBATCH --account=def-ben

module load nixpkgs/16.09  gcc/7.3.0
module load gmap-gsnap/2018-07-04
module load samtools/1.9

time gmap -D /home/songxy/scratch/tropicalis_transcriptome/tropicalis_genome/db_gmap_tropicalis_v91 -d db_gmap_tropicalis_v91 -A -B 5 -t 15 -f samse /home/songxy/scratch/tropicalis_transcriptome/transcriptome_building/tropicalis_transcriptome_trinityOut.Trinity.SuperTrans.fasta | samtools view -S -b > /home/songxy/scratch/tropicalis_transcriptome/mapping_transcriptome_to_genome/tropicalis_denovoT_tropicalisv91_genome_gmap.bam


```


