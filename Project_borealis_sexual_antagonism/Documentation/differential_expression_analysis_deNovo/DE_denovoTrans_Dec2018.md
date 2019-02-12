# Repeat DE for borealis transcriptome Dec2018
## transcriptome building
Graham has the newest version of Trinity, hence, I rebuild the transcriptome with the following bash script. 
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=1-22:00
#SBATCH --mem=210G
#SBATCH --job-name=borealis_transcriptome_dec2018
#SBATCH --account=def-ben

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
## Supertranscriptome
#### Perform DTU for supertranscriptome dec 2018
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --time=00:12:00
#SBATCH --mem=1G
#SBATCH --job-name=borealis_supertranscriptome_DTU
#SBATCH --account=def-ben


module load nixpkgs/16.09 intel/2016.4
#module load star/2.6.1a
module load hisat2/2.1.0
module load r/3.5.0
module load trinity/2.8.4
module load python/3.7.0


/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/intel2016.4/trinity/2.8.4/trinityrnaseq-Trinity-v2.8.4/Analysis/SuperTranscripts/DTU/dexseq_wrapper.pl --genes_fasta /home/songxy/projects/def-ben/songxy/borealis_transcriptome/data/borealis_denovo_transcriptome_dec2018_trinityOut/supertranscriptome_dec2018/borealis_transcriptome_trinityOut.Trinity.SuperTrans.fasta --genes_gtf /home/songxy/projects/def-ben/songxy/borealis_transcriptome/data/borealis_denovo_transcriptome_dec2018_trinityOut/supertranscriptome_dec2018/borealis_transcriptome_trinityOut.Trinity.SuperTrans.gtf --samples_file /home/songxy/projects/def-ben/songxy/borealis_transcriptome/data/trimmed_reads/trimmed_reads_individual/borealis_liver_trimmed_reads_samples.txt --CPU 10 --out_prefix /home/songxy/projects/def-ben/songxy/borealis_transcriptome/analysis/supertranscriptome_dec2018_DTU/DTU --aligner hisat2

```
## Transcriptome
##### Mapping borealis transcriptome to laevis genome
Tried it in graham
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=10:00:00
#SBATCH --mem=15G
#SBATCH --job-name=borealis_gmap_indexing
#SBATCH --account=def-ben

module load nixpkgs/16.09  gcc/7.3.0
module load gmap-gsnap/2018-07-04
module load samtools/1.9

time gmap -D /home/songxy/projects/def-ben/songxy/genome/laevis_genome/db_gmap_xl92 -d laevis92_gmap -A -B 5 -t 20 -f samse --cross-species /home/songxy/projects/def-ben/songxy/borealis_transcriptome/data/borealis_denovo_transcriptome_dec2018_trinityOut/borealis_transcriptome_trinityOut.Trinity.fasta | samtools view -S -b > /home/songxy/projects/def-ben/songxy/borealis_transcriptome/analysis/mapping_xb_denovoTrans_xl_genomev92_gmap/borealis_denovoT_laevisv92_genome_gmap.bam
```
didn't work because it said it didn't find the index file. Weird. Will try later when I have time. 

I did it in infor like this:
```bash
time gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92 -d laevis92_gmap -A -B 5 -t 25 -f samse --cross-species /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_dec2018/transcriptome/borealis_transcriptome_trinityOut.fasta | samtools view -S -b > /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_dec2018/analysis/transcriptome/mapping_xb_denovoTrans_xl_genomev92_gmap/borealis_denovoT_laevisv92_genome_gmap.bam
```
#### Collapsing the transcriptome 
Following the pipeline that I did in `Collapsing_transcripts.md`.

filter the mapping read remove the unmapped reads, which is indicated by flag 0x04. Then I piled the filtered output to bedtools, which extracts alignment coordinates for transcripts based on their CIGAR strings;
```
samtools view -F 0x04 -b borealis_denovoT_laevisv92_genome_gmap.bam | bedtools bamtobed -i > borealis_denovoT_laevisV92_genome_gmap_bedfile.bed
```
filter the bedfile and keep only the hit with the highest mapping quality score; if there are two hit with the same mapping quality score, keep the first one only;
```
perl ~/script/orthologs_identification_filterBedfile.pl borealis_denovoT_laevisV92_genome_gmap_bedfile.bed > borealis_de_laevisV92_genome_gmap_bedfile_filtered.tsv
```
##### Transcript expression levele Quantificaiton
Use Kallisto to quantify the expression count
```bash
#Indexing the transcriptomes:
kallisto index -i /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_dec2018/transcriptome/kallisto_index/borealis_transcriptome_trinityOut.fasta.kallisto_idx /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_dec2018/transcriptome/borealis_transcriptome_trinityOut.fasta

# Transcript abundance quantification:
kallisto quant -i /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta.kallisto_idx  -o female_rep7 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4082_girl_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4082_girl_liver_R2_scythe.fastq.gz)
# compile the read count for each samples into a matrix; used a script included in the Trinity RNAseq analysis package
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/util/abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix borealis_liver  --name_sample_by_basedir female_rep1/abundance.tsv female_rep2/abundance.tsv female_rep3/abundance.tsv female_rep4/abundance.tsv male_rep1/abundance.tsv male_rep2/abundance.tsv male_rep3/abundance.tsv male_rep4/abundance.tsv
```


