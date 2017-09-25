# Kallisto + EdgeR
Building abundance index and matrix with Kallisto with Trinity script. Trimmed reads were mapped to 
```
time perl /home/xue/trinityrnaseq-2.2.0/util/align_and_estimate_abundance.pl --transcripts /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta --seqType fa --right /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/borealis_allSamples_right.fastq.gz --left /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/borealis_allSamples_left.fastq.gz --est_method kallisto --output_dir /home/xue/borealis_DE --trinity_mode --prep_reference

time perl /home/xue/trinityrnaseq-2.2.0/util/abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix trans_counts  --name_sample_by_basedir BJE3897_mom/abundance.tsv BJE3929_boy/abundance.tsv BJE4009_girl/abundance.tsv BJE4017_boy/abundance.tsv BJE4039_boy/abundance.tsv BJE4072_girl/abundance.tsv BJE4082_girl/abundance.tsv dad/abundance.tsv

```
