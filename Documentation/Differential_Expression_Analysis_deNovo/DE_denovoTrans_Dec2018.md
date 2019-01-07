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

## Mapping borealis transcriptome to laevis genome
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

time gmap -D /home/songxy/projects/def-ben/songxy/genome/laevis_genome/db_gmap_xl92 -d laevis92_gmap -A -B 5 -t 20 -f samse --cross-species /home/songxy/projects/def-ben/songxy/borealis_transcriptome/borealis_denovo_transcriptome_dec2018_trinityOut/borealis_transcriptome_trinityOut.Trinity.fasta | samtools view -S -b > /home/songxy/projects/def-ben/songxy/borealis_transcriptome/mapping_xb_denovoTrans_xl_genomev92_gmap/borealis_denovoT_laevisv92_genome_gmap.bam
```
