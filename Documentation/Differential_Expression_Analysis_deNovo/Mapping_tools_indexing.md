# Mapping tools indexing
Here I will document the codes that I used for indexing reference genome/transcriptome for different alignment tools.

# GMAP - Indexing laevis genome (v91)  
Below is how I build laevis genome (v91) database on Graham for GMAP
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=20G
#SBATCH --job-name=tropicalis_gmap_indexing
#SBATCH --account=def-ben

module load nixpkgs/16.09  gcc/7.3.0
module load gmap-gsnap/2018-07-04

gmap_build -d db_gmap_tropicalis_v91/ -D /home/songxy/scratch/tropicalis_transcriptome/tropicalis_genome/db_gmap_tropicalis_v91 -g /home/songxy/scratch/tropicalis_transcriptome/tropicalis_genome/XT9_1.fa.gz
```
