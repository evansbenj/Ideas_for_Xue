## mapping
mapping super transcript to laveis genome
```
gmap -d laevis92_gmap/ -D /home/xue/genome_data/laevis_genome/db_gmap_xl92 -A -B 5 -t 24 -f samse --cross-species /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_supertranscriptApproach/SuperDuper.fasta | samtools view -S -b > /home/xue/borealis_transcriptome/borealis_supertranscriptome_2018/SuperDuper_gmap_laevisGenome.bam

```

## DE
Working directory:
```
/home/xue/borealis_transcriptome/borealis_supertranscriptome_2018/borealis_supertranscript_de
```
kallisto indexing
```
time perl /home/xue/software/trinityrnaseq-Trinity-v2.4.0/util/align_and_estimate_abundance.pl --transcripts /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_supertranscriptApproach/SuperDuper.fasta --seqType fa --samples_file /home/xue/borealis_transcriptome/borealis_supertranscriptome_2018/borealis_supertranscript_de/borealis_de_liver_samplefile.tsv --est_method kallisto --output_dir /home/xue/borealis_transcriptome/borealis_supertranscriptome_2018/borealis_supertranscript_de/kallisto_out/ --trinity_mode --prep_reference
```
convert to estimate matrix
```
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/util/abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix borealis_supertranscript_liver  --name_sample_by_basedir female_rep1/abundance.tsv female_rep2/abundance.tsv female_rep3/abundance.tsv female_rep4/abundance.tsv male_rep1/abundance.tsv male_rep2/abundance.tsv male_rep3/abundance.tsv male_rep4/abundance.tsv
```
DE using EdgeR
```
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix /home/xue/borealis_DE/liver_mvsf/kallisto_out/borealis_liver.counts.matrix --method edgeR --samples_file /home/xue/borealis_DE/liver_mvsf/edgeR_out/borealis_liver_samples_files.tsv
```
DE using DESeq2
```

```

