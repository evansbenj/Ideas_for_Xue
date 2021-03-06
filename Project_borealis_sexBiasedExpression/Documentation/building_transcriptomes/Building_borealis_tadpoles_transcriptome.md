# Processing of borealis tadpole RNAseq data
Here, I repeated the same step that I did for tropicalis. 

## RNAseq raw data
#### Downloading Raw RNAseq data from sequencing center
```
```
#### Upload the Raw RNAseq data to info
```
```
#### Do fastqc on raw data
We did it for all fastq.gz file in the folder
```
for i in *fastq.gz ; do fastqc --outdir=/home/xue/borealis_tadpoles_gonad_transcriptome_Feb2019/data/2019_XB_gonad_RNAseq/fastqc_raw $i ; done
```
and then download the file to a local folder and then inspect with a browser like this:

```
scp xue@info.mcmaster.ca:/home/xue/borealis_tadpoles_gonad_transcriptome_Feb2019/data/2019_XB_gonad_RNAseq/fastqc_raw/*html .
```


## Scythe - raw data
more detail about scythe: https://github.com/vsbuffalo/scythe. 

#### how to run Scythe on all the file with one command
```
for i in *fastq.gz ; do name=$(grep -o "XBO[0-9]*_R[0-9]" <(echo $i)); /home/xue/software/scythe-master/scythe -a /home/xue/software/scythe-master/illumina_adapters.fa -p 0.1 $i | gzip > /home/xue/borealis_tadpoles_gonad_transcriptome_Feb2019/data/scythed_data/$name\_scythe.fastq.gz; done
```
What are the parameters
- -a: adapter file in fasta format
- -p: the prior contamination rate is 0.05; BenF used 0.1.
- -o: output file name

Other options:
- -q: Illumina's quality scheme (pipeline > 1.3) is used. Sanger or Solexa (pipeline < 1.3) qualities can be specified
- -n: one can specify the minimum match length argument with -n <integer> and 
- -M: the minimum length of sequence (discarded less than or equal to this parameter) to keep after trimming with -M <integer>


#### run fastqc on scythe data
We did it for all fastq.gz file in the folder
```
for i in *fastq.gz ; do fastqc --outdir=/home/xue/borealis_tadpoles_gonad_transcriptome_Feb2019/data/scythed_data/fastqc_scythed/ $i ; done
```

and then download the file to a local folder and then inspect with a browser like this:

```
scp xue@info.mcmaster.ca:/home/xue/tropicalis_gonad_transcriptome_Dec2018/scythed_data/*html .
```

## Trimmomatic - 
run it like this if it is paired end
learn how to run trimmomatic and what the command means: http://www.usadellab.org/cms/?page=trimmomatic

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

#### running Trimmomatic on all scythed data
```bash
perl ~/script/run_trimmomatic.pl /home/xue/borealis_tadpoles_gonad_transcriptome_Feb2019/data/scythed_data/ /home/xue/borealis_tadpoles_gonad_transcriptome_Feb2019/data/trimmed_data/
```

#### run fastqc on Trimmomatic trimmed data
We did it for all fastq.gz file in the folder
```bash
for i in *fastq.gz ; do fastqc $i; done
```

## Trinity - Transcriptome assembly
More detail about trinity here: https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity

Before start building transcriptome, I will like to combine all the R1_paired together and all the R2_paired together by doing the below. The reason is that we want to build once single transcriptome that we can map read to instead of building individual ones and trying to find what is the equavalent transcripts in each transcriptome. 
```bash
cat *_R1_paired.fastq.gz > XBO_tad_R1.fastq.gz
cat *_R2_paired.fastq.gz > XBO_tad_R2.fastq.gz
```

I built the *de novo* transcriptome in Graham since Graham has the newest version of Trinity (2.8). Once it is build, I will transfer the transcriptome back to info for down-stream analysis. I moved the data to Graham by:
```
scp xue@info.mcmaster.ca:/home/xue/borealis_tadpoles_gonad_transcriptome_Feb2019/data/trimmed_data/*XBO_tad_* .
```

Graham alread have the newest version of Trinity (Trinity 2.8.4) installed, I am running the newest verion on Graham with the below bash script.
```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=5-12:00
#SBATCH --mem=230G
#SBATCH --job-name=borealis_transcriptome_Feb2019
#SBATCH --account=def-ben

module load nixpkgs/16.09  gcc/7.3.0 nixpkgs/16.09
module load openmpi/3.1.2
#module load trinity/2.8.4
module load samtools/1.9
module load salmon/0.11.3
module load bowtie2/2.3.4.3
module load jellyfish/2.2.6
module load trinity/2.8.4

Trinity --seqType fq --left /home/songxy/projects/def-ben/songxy/borealis_tad_gonad_transcriptome/data/trimmed_data/XBO_tad_R1.fastq.gz --right /home/songxy/projects/def-ben/songxy/borealis_tad_gonad_transcriptome/data/trimmed_data/XBO_tad_R2.fastq.gz --CPU 20 --full_cleanup --max_memory 200G --min_kmer_cov 2 --include_supertranscripts --output home/songxy/projects/def-ben/songxy/borealis_tad_gonad_transcriptome/data/borealis_tad_trinityOut
```
The one ran in info gave my error message and stopped running. The one in graham ran for a little bit more than 1 days with max 191Gb of memories. Much shorter time than I expected. 

#### output from trinity
- the output will be a fasta file by default. It will have the transcript name and transcript sequence for each transcript it builds.
- example of a transcript name is: >TRINITY_DN86389_c0_g1_i1 len=398 path=[0:0-159 2:160-208 3:209-275 5:276-397]
  - the c means cluster
  - the g means gene
  - the i means isoform: different isoform of a gene are supposely splice variants if Trinity did the assembly correctly.

## GMAP - Mapping to *X. laevis* genome (v92)
# Identify transcript's genomic location

## Mapping borealis transcriptome to laevis genome

Do mapping in graham since it is faster
```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --time=04:20:00
#SBATCH --mem=10G
#SBATCH --job-name=tropicalis_gmap_indexing
#SBATCH --account=def-ben

module load nixpkgs/16.09  gcc/7.3.0
module load gmap-gsnap/2018-07-04
module load samtools/1.9

time gmap -D /home/songxy/projects/def-ben/songxy/genome/laevis_genome/db_gmap_xl92 -d gmap_laevis_v92 -A -B 5 -t 15 -f samse /home/songxy/projects/def-ben/songxy/borealis_tad_gonad_transcriptome/data/transcriptome/borealis_tad_goand_transcriptome.fasta | samtools view -S -b > /home/songxy/projects/def-ben/songxy/borealis_tad_gonad_transcriptome/analysis/transcriptome/mapping_trans_laevisGenome92_gmap/borealisTad_denovoT_laevisGenome92_gmap.bam
```
Gmap on graham ran for 3.5hrs. 

## Coverting to bedfile
Following the pipeline that I did in `Collapsing_transcripts.md`.

I removed the unmapped reads, which is indicated by flag 0x04. Then I piped the filtered output to bedtools, which extracts alignment coordinates for transcripts based on their CIGAR strings;
```
#transcriptome
samtools view -F 0x04 -b borealisTad_denovoT_laevisGenome92_gmap.bam | samtools sort -n | bedtools bamtobed -i > borealisTad_denovoT_laevisGenome92_gmap_bedfile.bed

#supertranscriptome
samtools view -F 0x04 -b borealisTad_supertrans_laevisGenome92_gmap.bam | samtools sort -n | bedtools bamtobed -i > borealisTad_supertrans_laevisGenome92_gmap_bedfile.bed

```
With a perl script, I filterred the bedfile and keep only the hit with the highest mapping quality score; if there are two hit with the same mapping quality score, keep the first one only; **This might cause problem in the future**. 
  - one of the reason that the transcript mapped to multiple location with the same mapping quality is that the transcript was not assembled correctly. This transcript should be drop from the data set
  - another reason could be that it was mapping to the homologous chromosome. This can be fix by looking at the supertranscriptome and see which chromosome the corresponding supertranscript mapped to. -> **look at it in R** -> this can be used as a extra step to confirm genomic location in the case of doubt.
  
```
perl ~/script/orthologs_identification_filterBedfile.pl borealis_denovoT_laevisV92_genome_gmap_bedfile.bed > borealis_de_laevisV92_genome_gmap_bedfile_filtered.tsv
```
The path to bedfile is 
```
/home/xue/borealis_transcriptome/borealis_denovo_transcriptome_dec2018/analysis/transcriptome/mapping_xb_denovoTrans_xl_genomev92_gmap/
```

