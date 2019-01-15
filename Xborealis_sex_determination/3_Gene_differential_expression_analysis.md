# Gene differential expression_analysis
## Workflow 
- [ ] Collapsed transcriptome into gene level -> use supertranscriptome?
- [ ] Raw count quantafication - Kallisto
- [ ] Differential expression analysis - EdgeR

## Collapsing

## Kallisto
```
#Indexing the transcriptomes
kallisto index -i /home/xue/tropicalis_gonad_transcriptome_Dec2018/analysis/transcript_expression_raw_count/kallisto_indexing_transcriptome/tropicalis_transcriptome_trinityOut.Trinity.fasta.kallisto_idx /home/xue/tropicalis_gonad_transcriptome_Dec2018/data/tropicali_gonad_transcriptome_trinityOut/tropicalis_transcriptome_build_dec2018/tropicalis_transcriptome_trinityOut.Trinity.fasta

#Transcript abundance quantification - female
kallisto quant -i /home/xue/tropicalis_gonad_transcriptome_Dec2018/analysis/transcript_expression_raw_count/kallisto_indexing_transcriptome/tropicalis_transcriptome_trinityOut.Trinity.fasta.kallisto_idx  -o female_rep1 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3897_mom_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3897_mom_liver_R2_scythe.fastq.gz)

#Transcript abundance quantification - male
kallisto quant -i /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/mapped_borealis_transcripts/mapped_borealis_transcript_sequence.fa.kallisto_idx  -o male_rep4 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4039_boy_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4039_boy_liver_R2_scythe.fastq.gz)

#To compute the matrix (time cost: ~10min):
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/util/abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix borealis_liver  --name_sample_by_basedir female_rep1/abundance.tsv female_rep2/abundance.tsv female_rep3/abundance.tsv female_rep4/abundance.tsv male_rep1/abundance.tsv male_rep2/abundance.tsv male_rep3/abundance.tsv male_rep4/abundance.tsv
```
## EdgeR
