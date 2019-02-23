following tutorial: `https://github.com/trinityrnaseq/trinityrnaseq/wiki/DiffTranscriptUsage`

How to run it on info:
```
/home/xue/software/trinityrnaseq-Trinity-v2.7.0-PRERELEASE/Analysis/SuperTranscripts/DTU/dexseq_wrapper.pl \
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
#SBATCH --mem=100G
#SBATCH --job-name=borealis_supertranscriptome_dec2018
#SBATCH --account=def-ben

module load nixpkgs/16.09  gcc/7.3.0
module load nixpkgs/16.09  intel/2018.3
module load samtools/1.9
module load r/3.5.1
module load star/2.6.1c
module load python/3.7.0
module load trinity/2.8.4


/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/intel2016.4/trinity/2.8.4/trinityrnaseq-Trinity-v2.8.4/Analysis/SuperTranscripts/DTU/dexseq_wrapper.pl \
       --genes_fasta /home/songxy/projects/def-ben/songxy/borealis_transcriptome/supertranscriptome/trinity_genes.fasta \
       --genes_gtf /home/songxy/projects/def-ben/songxy/borealis_transcriptome/supertranscriptome/trinity_genes.gtf \
       --samples_file /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/borealis_liver_trimmed_reads_samples.txt \
       --CPU 20 \
       --out_prefix DTU --aligner STAR
```

```
#samples.txt

female female_rep1   /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE3897_mom_liver_R1_scythe.fastq.gz     /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE3897_mom_liver_R2_scythe.fastq.gz
female female_rep2   /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE4009_girl_liver_R1_scythe.fastq.gz     /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE4009_girl_liver_R2_scythe.fastq.gz
female female_rep3   /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE4072_girl_liver_R1_scythe.fastq.gz     /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE4072_girl_liver_R2_scythe.fastq.gz
female female_rep4   /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE4082_girl_liver_R1_scythe.fastq.gz     /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE4082_girl_liver_R2_scythe.fastq.gz
male male_rep1   /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE3896_dad_liver_R1_scythe.fastq.gz     /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE3896_dad_liver_R2_scythe.fastq.gz
male   male_rep2   /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE3929_boy_liver_R1_scythe.fastq.gz     /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE3929_boy_liver_R2_scythe.fastq.gz
male male_rep3   /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE4017_boy_liver_R1_scythe.fastq.gz     /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE4017_boy_liver_R2_scythe.fastq.gz
male male_rep4   /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE4039_boy_liver_R1_scythe.fastq.gz     /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE4039_boy_liver_R2_scythe.fastq.gz
```

## Salmon
#### indexing
```
salmon index -t /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_dec2018/supertranscriptome/borealis_superTrans.fasta -i /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_dec2018/analysis/supertranscriptome/count_salmon/borealis_superTrans_index

```
#### quantifying  
```bash
#!/bin/bash

index_dir=/home/xue/borealis_transcriptome/borealis_denovo_transcriptome_dec2018/analysis/supertranscriptome/count_salmon/borealis_superTrans_index

sample_dir=/home/xue/borealis_transcriptome/data/trimmed

out_dir=/home/xue/borealis_transcriptome/borealis_denovo_transcriptome_dec2018/analysis/supertranscriptome/count_salmon

cd /home/xue/borealis_transcriptome/data/trimmed

for i in *fastq.gz; do 
        sample_name=$(grep -o "BJE[0-9]*_[a-z]*" <(echo $i)) 
        
        salmon quant -i ${index_dir} -l A \
          -1 ${sample_dir}/${sample_name}_liver_R1_paired.fastq.gz \
          -2 ${sample_dir}/${sample_name}_liver_R2_paired.fastq.gz \
          -p 15 --validateMappings --rangeFactorizationBins 4 \
          --seqBias --gcBias \
          -o ${out_dir}/${sample_name}_quant
done
```



















