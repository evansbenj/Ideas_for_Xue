following tutorial: `https://github.com/trinityrnaseq/trinityrnaseq/wiki/DiffTranscriptUsage`

How to run it on info:
```
/home/xue/software/trinityrnaseq-Trinity-v2.8.4/Analysis/SuperTranscripts/DTU/dexseq_wrapper.pl \
       --genes_fasta /home/xue/borealis_transcriptome/borealis_supertranscriptome_2018/supertranscriptome_dec2018/trinity_genes.fasta \
       --genes_gtf /home/xue/borealis_transcriptome/borealis_supertranscriptome_2018/supertranscriptome_dec2018/trinity_genes.gtf \
       --samples_file /home/xue/borealis_transcriptome/borealis_supertranscriptome_2018/borealis_superTrans_DTU/borealis_de_liver_samplefile.tsv \
       --out_prefix DTU --aligner STAR
```

How to run it in Graham:
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

/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/intel2016.4/trinity/2.8.4/trinityrnaseq-Trinity-v2.8.4/Analysis/SuperTranscripts/$TRINITY_HOME/Analysis/SuperTranscripts/DTU/dexseq_wrapper.pl \
       --genes_fasta trinity_genes.fasta \
       --genes_gtf trinity_genes.gtf \
       --samples_file samples.txt \
       --out_prefix DTU --aligner STAR
```
