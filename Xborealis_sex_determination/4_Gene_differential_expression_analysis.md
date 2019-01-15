# Gene differential expression_analysis
## Workflow 
- [ ] Collapsed transcriptome into gene level -> use supertranscriptome?
- [ ] Raw count quantafication - Kallisto
- [ ] Differential expression analysis - EdgeR

## Collapsing

## Kallisto
```

time perl /home/xue/software/trinityrnaseq-Trinity-v2.4.0/util/align_and_estimate_abundance.pl --transcripts /home/xue/tropicalis_gonad_transcriptome_Dec2018/data/tropicali_gonad_transcriptome_trinityOut/tropicalis_transcriptome_build_dec2018/tropicalis_transcriptome_trinityOut.Trinity.fasta --seqType fa --samples_file /home/xue/borealis_DE/gonad_de_sex/borealis_gonad_kallisto_denovo/broealis_de_gonad_samplefile.tsv --est_method kallisto --output_dir /home/xue/borealis_DE/gonad_de_sex/borealis_gonad_kallisto_denovo/ --trinity_mode --prep_reference


#Indexing the transcriptomes
kallisto index -i /home/xue/tropicalis_gonad_transcriptome_Dec2018/analysis/transcript_expression_raw_count/kallisto_indexing_transcriptome/tropicalis_transcriptome_trinityOut.Trinity.fasta.kallisto_idx /home/xue/tropicalis_gonad_transcriptome_Dec2018/data/tropicali_gonad_transcriptome_trinityOut/tropicalis_transcriptome_build_dec2018/tropicalis_transcriptome_trinityOut.Trinity.fasta

#Transcript abundance quantification
for i in /home/xue/tropicalis_gonad_transcriptome_Dec2018/data/trim/*_R1_paired.fastq.gz; do name=$(grep -o "XT[0-9]*" <(echo $i));r1=/home/xue/tropicalis_gonad_transcriptome_Dec2018/data/trim/$name\_R1_paired.fastq.gz;r2=/home/xue/tropicalis_gonad_transcriptome_Dec2018/data/trim/$name\_R2_paired.fastq.gz; kallisto quant -i /home/xue/tropicalis_gonad_transcriptome_Dec2018/analysis/transcript_expression_raw_count/kallisto_indexing_transcriptome/tropicalis_transcriptome_trinityOut.Trinity.fasta.kallisto_idx  -o /home/xue/tropicalis_gonad_transcriptome_Dec2018/analysis/transcript_expression_raw_count/kallisto_indexing_transcriptome/$name <(gunzip -c $r1) <(gunzip -c $r2);done

#To compute the matrix (time cost: ~10min)-didnt do yet
#time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/util/abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix borealis_liver  --name_sample_by_basedir female_rep1/abundance.tsv female_rep2/abundance.tsv female_rep3/abundance.tsv female_rep4/abundance.tsv male_rep1/abundance.tsv male_rep2/abundance.tsv male_rep3/abundance.tsv male_rep4/abundance.tsv
```
## EdgeR
