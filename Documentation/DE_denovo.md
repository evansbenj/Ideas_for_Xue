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

Differential expression between two individuals using Trinity scripts. Specifying EdgeR as the method as it is the recommended method :
```
time perl /home/xue/trinityrnaseq-2.2.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix trans_counts.counts.matrix --method edgeR --dispersion 0.1

time perl /home/xue/trinityrnaseq-2.2.0/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ../trans_counts.TMM.EXPR.matrix -P 1e-3 -C 2
```
Output DE files were stored in: ``` /home/xue/borealis_DE/```.

# get the id of DE transcript (FDR<0.05)
```
awk '$5<0.05 {print $1}' all.counts.matrix.female_vs_male.edgeR.DE_results > all_mvsf_fdr005.tsv
awk '$2 > 2 && $5<0.05  {print $1}' gonad.counts.matrix.female_gonad_vs_male_gonad.edgeR.DE_results > gonad_fdr005.tsv
awk '$5<0.05 {print $1}' liver.counts.matrix.female_liver_vs_male_liver.edgeR.DE_results > liver_fdr005.tsv
awk '$5<0.05 {print $1}' tissue.counts.matrix.gonad_vs_liver.edgeR.DE_results > tissue_fdr005.tsv

```

# Extract sequences of DE genes
A script were used to extract the sequence of DE transcripts and output to a .fa file.
```
. get_trans_fdr005.sh
```
# Blast DE X.borealis sequence against annotated X.laevis genome
To find out the chromosomal location of the extracted X.borealis sequences, the extracted sequence will be mapped to the annotated X.laevis transcriptome using Blast. 
```
blastn -task blastn -db /home/xue/borealis_DE/xl_genome/db_Xlaevis_v91 -outfmt 6 -evalue 0.00005 -query  /home/xue/borealis_DE/all_mvsf/all_trans_fdr005.fa -out /home/xue/borealis_DE/all_mvsf/all_mvsf_blastout
```
I cleaned up the blastn output by filtering by e-value and bitscore.
```
perl filter_blastout.pl all_mvsf_blastout > tophit_ebs_blastout 
```
I extract the transcript ID of DE genes that is mapped to the Chr8L and is mapped to the first 50Mbp of chr8L.
```
awk '$2=="chr8L" && $9<50000000 && $10<50000000{print $1}' tophit_ebs_blastout > tophit_chr8L_50mb_blastout
```
The I extracted the DE information from EdgeR result for the above genes.
```
. get_DE_chr8L.sh
```

Then I examine the direction of gene expression by counting the number of DE have positive logFC and negative logFC. Positive logFC means the transcript has higher expression level in male and nagetive logFC means the transcript has higher expression level in female.
```
awk '$2>0{print}' DE_chr8L_50mb_EdgeR.tsv |wc -l
awk '$2<0{print}' DE_chr8L_50mb_EdgeR.tsv |wc -l
```
I combined all of the above into a bash script so that it will be easier to repeat the above in the future. I also trying to do the same in R, which seems like it will be easier (in progress).
```
. postDE_analysis.sh
```
The above steps are repeated for DE result from two comparision: male all sample vs female all sample; male liver vs female liver. The result are below:
- location distribution of DE genes


|chr #|all samples|liver samples|
|--- | --- | ---|
|total|957|432|
|chr8L|229 (23.9%)|160 (37%)|
|chr8S|70 (7.3%)|9%|
|chr2L|4.7%|3%|
|chr2S|5.5%|2.8%|

- DE genes that were mapped to chr8L in the first 50Mb have higher gene expression in which sex? positive logFC (>0) means the transcript has higher expression level in male and nagetive logFC (<0) means the transcript has higher expression level in female.


|logFC|all|liver
|--- | --- | ---|
|>0|34|24|
|<0|116 (77%)|92 (86%)|

From this point on, I will do the preceeding analysis with DE result from liver samples only to complete the pipline. 

# Binning trancripts 
- group overlapping transcripts into bins.  
- count the number of transcripts in each bin. If 100k transcripts were mapped to chr8L, how many of them overlapped?
- sum the expression level of transcripts in the same bin -> see if bins on chr8L non-recombinating region have different expression level than other bins (on chr8L non-sex-linked region and other autosomes)

# Compare to *X.laevis* chr8L
- do the same pipleline with XL RNA-seq data
- compare the expression level of DE genes in XB male, XB female, XL male, and XL female
- do a volcano plot to visualize those expression 

# Do trinity assembly for XB transcriptome with different setting






