# Collapsing transcripts using corset
Corset clusters transcripts from de novo transcriptome assembly to gene-level counts based on shared reads and expression patterns. 
- Corset documentation: https://github.com/Oshlack/Corset/wiki/InstallingRunningUsage
- Example: https://github.com/Oshlack/Corset/wiki/Example

### salmon
building index
```bash
salmon index --index Trinity --transcripts trinity_out_dir.Trinity.fasta
```
quantifying reads
```bash
##########bash script to run read quantificaitons for each individual separately using salmon  
#!/bin/bash

FILES=`ls /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE*_liver_R1_scythe.fastq.gz | sed 's/_liver_R1_scythe.fastq.gz//g'`
for F in $FILES ; do
        R1=${F}_liver_R1_scythe.fastq.gz
        R2=${F}_Rliver_R1_scythe.fastq.gz
        salmon quant --index Trinity --libType IU --dumpEq -1 $R1 -2 $R2 --output ${F}.out
done
```
### corset
install corset, which didn't go well on info. I will move to graham to do this analysis since Graham already have Corset 1.07 installed. 
```bash
./configure --prefix=/home/xue/software/Corset-version-1.07 --with-bam_inc=/usr/local/samtools/1.3.1/ --with-bam_lib=/usr/local/samtools/1.3.1/

make #make didn't work; 
```
run corset; in Graham, it ran for 50hrs, using max 7GB of memory. 
```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6 
#SBATCH --time=60:00:00
#SBATCH --mem=15G
#SBATCH --job-name=corset_grouping
#SBATCH --account=def-ben

module load nixpkgs/16.09  gcc/7.3.0
#module load corset/1.07
module load bowtie2/2.3.4.3
module load corset/1.07

corset -g 1,1,1,1,2,2,2,1 -n female_rep1,female_rep2,female_rep3,female_rep4,male_rep1,male_rep2,male_rep3,male_rep4 -f true -i salmon_eq_classes salmon_out/BJE*/aux_info/eq_classes.txt
```
Corset output a txt file, which contain the clustering of transcripts, and also a count table, which collapsed the expression count of transcripts that belong to the same cluster. The path to the files are:
```
/home/songxy/projects/def-ben/songxy/borealis_transcriptome/borealis_corset_collapsing/clusters.txt
/home/songxy/projects/def-ben/songxy/borealis_transcriptome/borealis_corset_collapsing/counts.txt
```


## Differential expression analysis






