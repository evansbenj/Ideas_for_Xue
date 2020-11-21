# Building borealis supertranscriptome

BenF already builded one supertranscriptome before but it is missing the gtf file, which I need as an input to the DTU script that trinity provided for DTU analysis. Hence, I decided to re-build this supertranscriptome following tutorial: `https://github.com/trinityrnaseq/trinityrnaseq/wiki/SuperTranscripts`.

The working directory is: /home/songxy/projects/def-ben/songxy/borealis_transcriptome/supertranscriptome/

I stated in info and try to generate Trinity SuperTranscripts like:
```bash
time /home/xue/software/trinityrnaseq-Trinity-v2.8.4/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py --trinity_fasta /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/trinity_out_dir.Trinity.fasta 
```
This script didn't run in info and always return error about unable to import python module. I decided to try to run it in Graham and below is the bash script. Since the Python script for building super-transcriptome is in the Trinity installation folder, we need to find the folder path like this:
```
module show Trinity
```
The path to the Trinity build is `/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/intel2016.4/trinity/2.8.4/trinityrnaseq-Trinity-v2.8.4`
Then we used a bash script to run the python script that build the super transcriptome
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3 
#SBATCH --time=01:12:00
#SBATCH --mem=10G
#SBATCH --job-name=borealis_supertranscriptome_dec2018
#SBATCH --account=def-ben

module load nixpkgs/16.09  gcc/7.3.0
#module load gnureadline-6.3.8 pip-18.0 setuptools-40.0.0 virtualenv-16.0.0 wheel-0.31.1
#module load trinity/2.8.4
#module load samtools/1.9
module load trinity/2.8.4
module load python/3.7.0

/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/intel2016.4/trinity/2.8.4/trinityrnaseq-Trinity-v2.8.4/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py --trinity_fasta trinity_out_dir.Trinity.fasta
~                   

```
and this generates two output files:
```
trinity_genes.fasta   :supertranscripts in fasta format
trinity_genes.gtf     :transcript structure annotation in gtf format
```

## Supertranscriptome - Dec 2018
I was rebuilding the borealis transcriptome with the newest version of Trinity on Dec 2018. The newest version included a flag `--include_supertranscripts` and I included this flag to get a supertranscriptome with the transcriptome. Since it was hard to install Trinity version 2.8.4 on info, I did the transriptome assembly in Graham. 
```bash
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
The path to the supertranscriptome is:
- on graham: 
  ```
  #sequence fasta
  /home/songxy/projects/def-ben/songxy/borealis_transcriptome/data/supertranscriptome_dec2018/borealis_superTrans.fasta
  
  #gene_trans_map
  /home/songxy/projects/def-ben/songxy/borealis_transcriptome/data/supertranscriptome_dec2018/borealis_superTrans_fasta.gene_trans_map
  
  #gtf
  /home/songxy/projects/def-ben/songxy/borealis_transcriptome/data/supertranscriptome_dec2018/borealis_superTrans.gtf
  ```

- on info: 
  ```
  #sequence fasta
  /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_dec2018/denovo_transcriptome/supertranscriptome/borealis_superTrans.fasta
  
  #gene_trans_map
  /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_dec2018/denovo_transcriptome/borealis_superTrans_fasta.gene_trans_map
  
  #gtf
  /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_dec2018/denovo_transcriptome/borealis_superTrans.gtf
  
  ```


