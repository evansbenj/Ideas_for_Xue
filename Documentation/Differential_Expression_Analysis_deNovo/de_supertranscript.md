## mapping
mapping super transcript to laveis genome
```
gmap -d laevis92_gmap/ -D /home/xue/genome_data/laevis_genome/db_gmap_xl92 -A -B 5 -t 24 -f samse --cross-species /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_supertranscriptApproach/SuperDuper.fasta | samtools view -S -b > /home/xue/borealis_transcriptome/borealis_supertranscriptome_2018/SuperDuper_gmap_laevisGenome.bam

```
convert to bedfile
```
samtools view -F 0x04 -b SuperDuper_gmap_laevisGenome.bam | bedtools bamtobed -i > SuperDuper_gmap_laevisGenome.bed
```
filter bedfile
```
perl ~/script/orthologs_identification_filterBedfile.pl SuperDuper_gmap_laevisGenome.bed > SuperDuper_gmap_laevisGenome_bedfile_filtered.tsv
```

## DE
Working directory:
```
/home/xue/borealis_transcriptome/borealis_supertranscriptome_2018/borealis_supertranscript_de
```
Kallisto
```
/home/unix/bhaas/GITHUB/trinityrnaseq/util/align_and_estimate_abundance.pl --transcripts Trinity.fasta --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference

build index
 time perl /home/xue/software/trinityrnaseq-Trinity-v2.4.0/util/align_and_estimate_abundance.pl --transcripts /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_supertranscriptApproach/SuperDuper.fasta   --est_method kallisto --output_dir /home/xue/borealis_transcriptome/borealis_supertranscriptome_2018/borealis_supertranscript_de/kallisto_out/ --trinity_mode --prep_reference
 
 kallisto
 time perl /home/xue/software/trinityrnaseq-Trinity-v2.4.0/util/align_and_estimate_abundance.pl --transcripts /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_supertranscriptApproach/SuperDuper.fasta --seqType fa --samples_file /home/xue/borealis_transcriptome/borealis_supertranscriptome_2018/borealis_supertranscript_de/borealis_de_liver_samplefile.tsv --est_method kallisto --output_dir /home/xue/borealis_transcriptome/borealis_supertranscriptome_2018/borealis_supertranscript_de/kallisto_out/ --trinity_mode
 
 time perl /home/xue/software/trinityrnaseq-2.2.0/util/align_and_estimate_abundance.pl --transcripts /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_supertranscriptApproach/SuperDuper.fasta --seqType fa --samples_file /home/xue/borealis_transcriptome/borealis_supertranscriptome_2018/borealis_supertranscript_de/borealis_de_liver_samplefile.tsv --est_method kallisto --output_dir /home/xue/borealis_transcriptome/borealis_supertranscriptome_2018/borealis_supertranscript_de/kallisto_out/ --trinity_mode
 
 #running kallisto manually
 kallisto quant -i /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_supertranscriptApproach/SuperDuper.fasta.kallisto_idx  -o female_rep1 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3897_mom_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3897_mom_liver_R2_scythe.fastq.gz)

kallisto quant -i /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_supertranscriptApproach/SuperDuper.fasta.kallisto_idx  -o male_rep4 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4039_boy_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4039_boy_liver_R2_scythe.fastq.gz)

female  female_rep1     /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3897_mom_liver_R1_scythe.fastq.gz   
female  female_rep2     /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4009_girl_liver_R1_scythe.fastq.gz 
female  female_rep3     /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4072_girl_liver_R1_scythe.fastq.gz 
female  female_rep4     /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4082_girl_liver_R1_scythe.fastq.gz 
male    male_rep1       /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3896_dad_liver_R1_scythe.fastq.gz   
male    male_rep2       /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3929_boy_liver_R1_scythe.fastq.gz   
male    male_rep3       /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4017_boy_liver_R1_scythe.fastq.gz   
male    male_rep4       /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4039_boy_liver_R1_scythe.fastq.gz   
```
convert to estimate matrix
```
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/util/abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix borealis_supertranscript_liver  --name_sample_by_basedir female_rep1/abundance.tsv female_rep2/abundance.tsv female_rep3/abundance.tsv female_rep4/abundance.tsv male_rep1/abundance.tsv male_rep2/abundance.tsv male_rep3/abundance.tsv male_rep4/abundance.tsv
```
DE using EdgeR
```
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix /home/xue/borealis_transcriptome/borealis_supertranscriptome_2018/borealis_supertranscript_de/borealis_supertranscript_liver.counts.matrix --method edgeR --samples_file /home/xue/borealis_transcriptome/borealis_supertranscriptome_2018/borealis_supertranscript_de/borealis_de_liver_samplefile.tsv
```
DE using DESeq2
```
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix /home/xue/borealis_transcriptome/borealis_supertranscriptome_2018/borealis_supertranscript_de/borealis_supertranscript_liver.counts.matrix --method DESeq2 --samples_file /home/xue/borealis_transcriptome/borealis_supertranscriptome_2018/borealis_supertranscript_de/borealis_de_liver_samplefile.tsv

```


## Summarize DE output
EdgeR result
```
/home/xue/borealis_transcriptome/borealis_supertranscriptome_2018/borealis_supertranscript_de/
```
Deseq2 result
```
awk ' $7<0.05  {print }' borealis_supertranscript_liver.counts.matrix.female_vs_male.DESeq2.DE_results > supertranscript_deseq2_fdr005.tsv

awk '($4 < -2||$4 >2) && $7<0.05  {print }' borealis_supertranscript_liver.counts.matrix.female_vs_male.DESeq2.DE_results > supertranscript_logfc2_fdr005.tsv
```

# Differential Transcript Usage
Working directory `/home/xue/borealis_transcriptome/borealis_supertranscriptome_2018/borealis_superTrans_DTU`. 

Following tutorial in https://github.com/trinityrnaseq/trinityrnaseq/wiki/DiffTranscriptUsage#differential-transcript-usage-via-supertranscripts. Firstly, I need to convert gff file to gtf file. 
```bash
#covert gff to gtf
gffread SuperDuper.gff -T -o SuperDuper.gtf
```
The gtf file is not complete and have the entire transcript as one exon. Need to re-run the supertranscriptome assembly to get a better gff file. 
```
```

```bash
~/software/trinityrnaseq-Trinity-v2.7.0-PRERELEASE/Analysis/SuperTranscripts/DTU/dexseq_wrapper.pl \
       --genes_fasta /home/xue/borealis_transcriptome/borealis_supertranscriptome_2018/SuperDuper.fasta \
       --genes_gtf /home/xue/borealis_transcriptome/borealis_supertranscriptome_2018/SuperDuper.gff \
       --samples_file /home/xue/borealis_transcriptome/borealis_supertranscriptome_2018/borealis_superTrans_DTU/borealis_de_liver_samplefile.tsv \
       --out_prefix borealis_DTU --aligner STAR
```

