# Differential Expression Analysis with Kallisto + EdgeR
Building abundance index and matrix with Kallisto using Trinity scripts. Each pair of trimmed reads (R1+R2) were mapped to the assembled transcriptome with reads from all the samples
```
time perl /home/xue/software/trinityrnaseq-Trinity-v2.5.0/util/align_and_estimate_abundance.pl --transcripts /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta --seqType fa --samples_file all_mvsf_samplefile.txt --est_method kallisto --output_dir /home/xue/borealis_DE/all_mvsf --trinity_mode --prep_reference


time perl /home/xue/software/trinityrnaseq-Trinity-v2.5.0/util/align_and_estimate_abundance.pl --transcripts /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta --seqType fa --samples_file gonad_mvsf_samplefile.tsv --est_method kallisto --output_dir /home/xue/borealis_DE/gonad_mvsf --trinity_mode --prep_reference

time perl /home/xue/software/trinityrnaseq-Trinity-v2.5.0/util/align_and_estimate_abundance.pl --transcripts /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta --seqType fa --samples_file liver_mvsf_samplefile.tsv --est_method kallisto --output_dir /home/xue/borealis_DE/liver_mvsf --trinity_mode --prep_reference


time perl /home/xue/software/trinityrnaseq-Trinity-v2.5.0/util/align_and_estimate_abundance.pl --transcripts /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta --seqType fa --samples_file liver_vs_gonad_samplefile.tsv --est_method kallisto --output_dir /home/xue/borealis_DE/liver_vs_gonad --trinity_mode --prep_reference

```
To compute the matrix (time cost: ~10min):
```
time perl /home/xue/trinityrnaseq-2.2.0/util/abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix all  --name_sample_by_basedir male_rep1/abundance.tsv male_rep2/abundance.tsv male_rep3/abundance.tsv male_rep4/abundance.tsv male_rep5/abundance.tsv male_rep6/abundance.tsv male_rep7/abundance.tsv female_rep1/abundance.tsv female_rep2/abundance.tsv female_rep3/abundance.tsv female_rep4/abundance.tsv female_rep5/abundance.tsv female_rep6/abundance.tsv female_rep7/abundance.tsv

time perl /home/xue/trinityrnaseq-2.2.0/util/abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix trans_counts  --name_sample_by_basedir BJE3896_dad_liver/abundance.tsv BJE3896_dad_testis/abundance.tsv BJE3897_mom_liver/abundance.tsv BJE3929_boy_liver/abundance.tsv BJE3929_boy_testis/abundance.tsv BJE4009_girl_liver/abundance.tsv BJE4009_girl_oviduct/abundance.tsv BJE4017_boy_liver/abundance.tsv BJE4017_boy_testis/abundance.tsv BJE4039_boy_liver/abundance.tsv BJE4039_boy_testis/abundance.tsv BJE4072_girl_liver/abundance.tsv BJE4072_girl_oviduct/abundance.tsv BJE4082_girl_liver/abundance.tsv BJE4082_girl_oviduct/abundance.tsv
```

Each takes about 45min - 1hr. The trinity perl script would CMD kallisto as following:
```
kallisto quant -i /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta.kallisto_idx  -o female_rep7 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4082_girl_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4082_girl_liver_R2_scythe.fastq.gz)
```

Differential expression between two individuals using Trinity scripts. Specifying EdgeR as the method as it is the recommended method when no biological replicate was available:
```
time perl /home/xue/trinityrnaseq-2.2.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix trans_counts.counts.matrix --method edgeR --dispersion 0.1

time perl /home/xue/trinityrnaseq-2.2.0/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ../trans_counts.TMM.EXPR.matrix -P 1e-3 -C 2
```
Output DE files were stored in: ``` /home/xue/borealis_DE/```.

# Extract sequences of DE genes
A script were used to extract the sequence of DE transcripts and output to a .fa file.
```
```
# Blast extracted sequence against annotated X.laevis transcriptome
To find out the chromosomal location of the extracted X.borealis sequences, the extracted sequence will be mapped to the annotated X.laevis transcriptome using Blast. 
```
```


