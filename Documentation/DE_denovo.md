# Differential Expression Analysis with Kallisto + EdgeR
Building abundance index and matrix with Kallisto using Trinity scripts. Each pair of trimmed reads (R1+R2) were mapped to the assembled transcriptome with reads from all the samples
```
time perl /home/xue/software/trinityrnaseq-2.2.0/util/align_and_estimate_abundance.pl --transcripts /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta --seqType fa --samples_file all_mvsf_samplefile.txt --est_method kallisto --output_dir /home/xue/borealis_DE/all_mvsf --trinity_mode --prep_reference


time perl /home/xue/software/trinityrnaseq-2.2.0/util/align_and_estimate_abundance.pl --transcripts /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta --seqType fa --samples_file gonad_mvsf_samplefile.tsv --est_method kallisto --output_dir /home/xue/borealis_DE/gonad_mvsf --trinity_mode --prep_reference

time perl /home/xue/software/trinityrnaseq-2.2.0/util/align_and_estimate_abundance.pl --transcripts /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta --seqType fa --samples_file liver_mvsf_samplefile.tsv --est_method kallisto --output_dir /home/xue/borealis_DE/liver_mvsf --trinity_mode --prep_reference


time perl /home/xue/software/trinityrnaseq-2.2.0/util/align_and_estimate_abundance.pl --transcripts /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta --seqType fa --samples_file liver_vs_gonad_samplefile.tsv --est_method kallisto --output_dir /home/xue/borealis_DE/liver_vs_gonad --trinity_mode --prep_reference

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
time perl /home/xue/software/trinityrnaseq-2.2.0//Analysis/DifferentialExpression/run_DE_analysis.pl --matrix liver.counts.matrix --method edgeR --samples_file liver_mvsf_samplefile_EdgeR.tsv

time perl /home/xue/trinityrnaseq-2.2.0/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ../trans_counts.TMM.EXPR.matrix -P 1e-3 -C 2
______________________
time perl /home/xue/software/trinityrnaseq-Trinity-v2.5.1//Analysis/DifferentialExpression/run_DE_analysis.pl --matrix ../liver.counts.matrix --method edgeR --samples_file ../liver_mvsf_samplefile_EdgeR.tsv
```


Output DE files were stored in: ``` /home/xue/borealis_DE/```.

# Get the id of DE transcript (FDR<0.05)
```
awk '$5<0.05 {print $1}' all.counts.matrix.female_vs_male.edgeR.DE_results > all_mvsf_fdr005.tsv
awk '$2 > 1 && $5<0.05  {print $1}' gonad.counts.matrix.female_gonad_vs_male_gonad.edgeR.DE_results > gonad_fdr005.tsv
awk '($2 < -1||$2 >1) && $5<0.05  {print $1}' liver.counts.matrix.female_liver_vs_male_liver.edgeR.DE_results > liver_fdr005.tsv
awk '$5<0.05 {print $1}' tissue.counts.matrix.gonad_vs_liver.edgeR.DE_results > tissue_fdr005.tsv

```

# Extract sequences of DE genes
A script were used to extract the sequence of DE transcripts and output to a fasta (.fa) file.
```
. get_trans_fdr005.sh
```
# Mapping Method 1: Blast DE X.borealis sequence against annotated X.laevis genome
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
Check the following to make sure the data are not biased due to the present of missing data.
- missing information, ex, 0 experssion level in male but have expression leve in female
- check how many DE mapped to chr8L or chr8S with same e-value (but different bitscore)

# Mapping with Star

Mapping assembled X.borealis transcriptome to X.laevis genome using STAR. First step is indexing. I used the default setting and provided gft file of the X.laevis genome (Coverting GFF3 file into GTF file):
```
/home/xue/borealis_DE/ref_approach/gffread/gffread/gffread -E XENLA_9.2_Xenbase.gff3 -T -o XENLA_9.2_Xenbase.gtf 

STAR --runThreadN 20 --runMode genomeGenerate --genomeDir /home/xue/borealis_DE/xl_genome/STAR_XLgenome_wGTF --genomeFastaFiles ../XL9_2.fa --sjdbGTFfile /home/xue/borealis_DE/ref_approach/XENLA_9.2_Xenbase.gtf --genomeChrBinNbits 16
```
Indexing without gft file:
```
STAR --runThreadN 20 --runMode genomeGenerate --genomeDir /home/xue/borealis_DE/xl_genome/STAR_XLgenome --genomeFastaFiles XL9_2.fa --genomeChrBinNbits 16
```
Second step is mapping.Mapping liver DE transcripts to XL genome:
```
STAR --runThreadN 20 --genomeDir /home/xue/borealis_DE/xl_genome/STAR_XLgenome --readFilesIn ../liver_trans_fdr005.fa  --genomeLoad LoadAndKeep --outFileNamePrefix ./Liver_StarOut --outSAMtype BAM Unsorted
```
Then I mapped the assembled X.borealis transcriptome to the X.laevis genome. Path to mapping output: ```/home/xue/borealis_DE/xbxl_mapping_STAR/```
```
time STAR --runThreadN 20 --genomeDir /home/xue/borealis_DE/xl_genome/STAR_XLgenome --readFilesIn /home/xue/borealis_DE/xl_genome/xb_transcriptome_subset/xb_tra_subset0.fa --genomeLoad LoadAndRemove -outFileNamePrefix /home/xue/borealis_DE/xl_genome/STAR_XLgenome/xb_STAR_1 --outSAMtype BAM Unsorted

time STAR --runThreadN 20 --genomeDir /home/xue/borealis_DE/xl_genome/STAR_XLgenome --readFilesIn /home/xue/borealis_DE/xl_genome/xb_transcriptome_subset/xb_tra_subset0.fa  --genomeLoad LoadAndRemove -outFileNamePrefix ./xb_STAR4 --outSAMtype BAM Unsorted

STAR --runThreadN 20 --genomeDir /home/xue/borealis_DE/xl_genome/STAR_XLgenome --readFilesIn /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta  --genomeLoad LoadAndKeep --outFileNamePrefix ./xbxl_ --outSAMtype BAM Unsorted

STAR --runThreadN 20 --genomeDir /home/xue/borealis_DE/xl_genome/STAR_XLgenome --readFilesIn ../xb_transcriptome_trinityout.fasta  --genomeLoad LoadAndKeep --outFileNamePrefix ./xbxl_ --outSAMtype BAM Unsorted
```


# Binning trancripts 
The Perl script that do the binning would: 
  - group overlapping transcripts into bins.  
  - count the number of transcripts in each bin. If 100k transcripts were mapped to chr8L, how many of them overlapped?
  - sum the expression level of transcripts in the same bin -> see if bins on chr8L non-recombinating region have different expression level than other bins (on chr8L non-sex-linked region and other autosomes)
  


# Compare to *X.laevis* chr8L
- do the same pipleline with XL RNA-seq data
- compare the expression level of DE genes in XB male, XB female, XL male, and XL female
- do a volcano plot to visualize those expression levels 

# Do trinity assembly for XB transcriptome with different setting
how did others do polyploid transcriptome assembly using Trinity? -> mostly default setting; 
```
Duan, J., Xia, C., Zhao, G., Jia, J. &Kong, X. Optimizing de novo common wheat transcriptome assembly using short-read RNA-Seq data. BMC Genomics 13, 392 (2012).
> ALLPATHSLG error correction, and the paired fragment length was set to 200 bp

Chopra, R. et al. Comparisons of De Novo Transcriptome Assemblers in Diploid and Polyploid Species Using Peanut (Arachis spp.) RNA-Seq Data.
> the default kmer of 25, minimum coverage of 2

Nakasugi, K., Crowhurst, R., Bally, J., Waterhouse, P. &Mittapalli, O. Combining Transcriptome Assemblies from Multiple De Novo Assemblers in the Allo-Tetraploid Plant Nicotiana benthamiana. PLoS One 9, (2014).
> default setting

Gutierrez-Gonzalez, J. J., Tu, Z. J., & Garvin, D. F. (2013). Analysis and annotation of the hexaploid oat seed transcriptome. BMC genomics, 14(1), 471.
> minimum assembled transcript length = 100 nt

Nakasugi, K., Crowhurst, R., Bally, J., Waterhouse, P. &Mittapalli, O. Combining Transcriptome Assemblies from Multiple De Novo Assemblers in the Allo-Tetraploid Plant Nicotiana benthamiana. PLoS One 9, (2014).
> Reads were normalized by the perl script insilico_read_normalization.pl (Trinity package) for Trinity;
> Two settings of k-mer 25 and k-mer 31

He, B. et al. Optimal assembly strategies of transcriptome related to ploidies of eukaryotic organisms. (2011). doi:10.1186/s12864-014-1192-7
> Trinity_release_20131110 (http://trinityrnaseq.github.io/) [16], which used in default parameter, kmer =25, was used in SASP strategy and its Command-line parameters were “--seqType fq --left Reads_1.fq --right Reads_2.fq --CPU 20”

McKain, M. R. et al. A phylogenomic assessment of ancient polyploidy and genome evolution across the Poales. Genome Biol. Evol. 8, evw060 (2016).
> Trinity (Release 2012-06-08; Grabherr et al. 2011 ) with default parameters

Li, H.-Z. et al. Evaluation of Assembly Strategies Using RNA-Seq Data Associated with Grain Development of Wheat (Triticum aestivum L.). PLoS One 8, e83530 (2013).
> Trinity (release 20120608) with minimum contig length  =  100, default parameter

Chen, S., McElroy, J. S., Dane, F. &Peatman, E. Optimizing Transcriptome Assemblies for Leaf and Seedling by Combining Multiple Assemblies from Three De Novo Assemblers. Plant Genome 8, 0 (2015).
> normalized with Trinity’s in silico read normalization (http://Trinityrnaseq.sourceforge.net/), with maximum coverage of 30 and k-mer size of 25
>  The merged assembly was first processed by CD-HIT-EST v4.6.1-2012-08-27 (http://weizhong-lab.ucsd.edu/cd-hit/) with the strictest parameters “-c 1.0-n 10” to remove identical fragments and then subjected to the EvidentialGene tr2aacds pipeline (http://arthropods.eugenes.org/EvidentialGene/about/EvidentialGene_trassembly_pipe.html). 

```
Try #1: increase min contig length to 300 -> stopped and abandoned
- By increasing the filter limit for minimum contig length (--min_contig_length), short fragment that is smaller than the limit would not be output by Trinity. Default is 200 so I increase it to 300. However, this might lead to some information loss.
```
/home/xue/software/trinityrnaseq-Trinity-v2.5.1/Trinity --seqType fq --left /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/*R1*  --right /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/*R2*  --CPU 20 --full_cleanup --max_memory 200G --min_contig_length 300 
```

Try #2: change the kmer size to 31 -> still in progress
- Maximum kmer size for Trinity is 32. One of the above paper set it to 31. 
Accouding to this paper (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5237644/), assembly produces less transcripts with larger kmer size. 
```
/home/xue/software/trinityrnaseq-Trinity-v2.5.1/Trinity --seqType fq --left /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/*R1*  --right /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/*R2*  --CPU 20 --full_cleanup --max_memory 200G --min_kmer_cov 31 --output /home/xue/xborealis_alltogether/xb_kmer31_trinity
```

Try #3: with RNA-seq from liver only -> stopped and abandoned 
- Since sex-specific genes could have different expression level in liver and gonad. When doing expression analysis, 
```
/home/xue/software/trinityrnaseq-Trinity-v2.5.1/Trinity --seqType fq --left /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/*liver_R1*  --right /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/*liver_R2*  --CPU 20 --full_cleanup --max_memory 200G  --output /home/xue/xborealis_alltogether/xb_liver_trinity
```
# count the number of zero expresssion combinations
It was done with a perl script.
- all

| |Female-0|Female-1|Female-2|Female-3|Female-4|
|---|---|---|---|---|---|
|Male-0 |96655  |52399  |17786  |7466   |2211|
|Male-1 |33485  |41692  |38417  |29014  |14119|
|Male-2 |11046  |33589  |55878  |69523  |56068|
|Male-3 |4474   |22143  |60984  |125930 |167275|


- DE 

| |Female-0|Female-1|Female-2|Female-3|Female-4|
|---|---|---|---|---|---|
|Male-0|17|2|5|7|65|
|Male-1|16|0|0|3|25|
|Male-2|13|0|0|0|0|
|Male-3|35|2|0|0|0|

- on chr8L


- on chr8L sex link region

 


# things that we could do but not urgent
- repeat DE analysis with different softwware package (RESM, DESeq2 etc)
- repeat the pipeline with DE comparison pair (all male smaples vs all female samples, male gonad vs female gonad, liver vs gonad)

# Path to laevis RNAseq read data
```
/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Session-SRA
```









