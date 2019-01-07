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
## Perform DTU for supertranscriptome dec 2018
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


/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/intel2016.4/trinity/2.8.4/trinityrnaseq-Trinity-v2.8.4/Analysis/SuperTranscripts/DTU/dexseq_wrapper.pl \
               --genes_fasta /home/songxy/projects/def-ben/songxy/borealis_transcriptome/supertranscriptome/trinity_genes.fasta \
                      --genes_gtf /home/songxy/projects/def-ben/songxy/borealis_transcriptome/supertranscriptome/trinity_genes.gtf \
                             --samples_file /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/borealis_liver_trimmed_reads_samples.txt \
                                    --CPU 10 \
                                           --out_prefix DTU --aligner hisat2

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
