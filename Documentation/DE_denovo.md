# Differential Expression Analysis with Kallisto + EdgeR
Building transcriptome index and abundance matrix with Kallisto using Trinity perl scripts. Each pair of trimmed reads (R1+R2) were pseudoaligned to the assembled transcriptome that were built with RNA-seq reads from all X.borealis individuals. More about how kallisto works, please see my note here (https://github.com/Sarahxys/Borealis/blob/master/About%20differential%20expression%20analysis.md). 
```
time perl /home/xue/software/trinityrnaseq-2.2.0/util/align_and_estimate_abundance.pl --transcripts /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta --seqType fa --samples_file all_mvsf_samplefile.txt --est_method kallisto --output_dir /home/xue/borealis_DE/all_mvsf --trinity_mode --prep_reference


time perl /home/xue/software/trinityrnaseq-2.2.0/util/align_and_estimate_abundance.pl --transcripts /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta --seqType fa --samples_file gonad_mvsf_samplefile.tsv --est_method kallisto --output_dir /home/xue/borealis_DE/gonad_mvsf --trinity_mode --prep_reference

time perl /home/xue/software/trinityrnaseq-2.2.0/util/align_and_estimate_abundance.pl --transcripts /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta --seqType fa --samples_file liver_mvsf_samplefile.tsv --est_method kallisto --output_dir /home/xue/borealis_DE/liver_mvsf --trinity_mode --prep_reference


time perl /home/xue/software/trinityrnaseq-2.2.0/util/align_and_estimate_abundance.pl --transcripts /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta --seqType fa --samples_file liver_vs_gonad_samplefile.tsv --est_method kallisto --output_dir /home/xue/borealis_DE/liver_vs_gonad --trinity_mode --prep_reference

```
The Trinity perl script CMD-ed kallisto as following. Kallisto can do quantification for one sample at a time. Each quantification takes about 45min - 1hr. I would run the same command if I were to run kallisto by myself. One extra thing that I can do if I am running the quantification step without the Trinity scripts is that I can do bootstraping with 
```
#Indexing the transcriptomes
kallisto index -i /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta.kallisto_idx /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta

#Transcript abundance quantification
kallisto quant -i /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta.kallisto_idx  -o female_rep7 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4082_girl_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4082_girl_liver_R2_scythe.fastq.gz)
```

To compute the matrix (time cost: ~10min):
```
time perl /home/xue/trinityrnaseq-2.2.0/util/abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix all  --name_sample_by_basedir male_rep1/abundance.tsv male_rep2/abundance.tsv male_rep3/abundance.tsv male_rep4/abundance.tsv male_rep5/abundance.tsv male_rep6/abundance.tsv male_rep7/abundance.tsv female_rep1/abundance.tsv female_rep2/abundance.tsv female_rep3/abundance.tsv female_rep4/abundance.tsv female_rep5/abundance.tsv female_rep6/abundance.tsv female_rep7/abundance.tsv

time perl /home/xue/trinityrnaseq-2.2.0/util/abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix trans_counts  --name_sample_by_basedir BJE3896_dad_liver/abundance.tsv BJE3896_dad_testis/abundance.tsv BJE3897_mom_liver/abundance.tsv BJE3929_boy_liver/abundance.tsv BJE3929_boy_testis/abundance.tsv BJE4009_girl_liver/abundance.tsv BJE4009_girl_oviduct/abundance.tsv BJE4017_boy_liver/abundance.tsv BJE4017_boy_testis/abundance.tsv BJE4039_boy_liver/abundance.tsv BJE4039_boy_testis/abundance.tsv BJE4072_girl_liver/abundance.tsv BJE4072_girl_oviduct/abundance.tsv BJE4082_girl_liver/abundance.tsv BJE4082_girl_oviduct/abundance.tsv
```



Differential expression between two conditions (ex, male_liver samples vs female_liver samples) with EdgeR using Trinity scripts. Specifying EdgeR as the method as it is the recommended method :
```
time perl /home/xue/software/trinityrnaseq-2.2.0//Analysis/DifferentialExpression/run_DE_analysis.pl --matrix liver.counts.matrix --method edgeR --samples_file liver_mvsf_samplefile_EdgeR.tsv

time perl /home/xue/trinityrnaseq-2.2.0/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ../trans_counts.TMM.EXPR.matrix -P 1e-3 -C 2

#with a newer version of Trinity
time perl /home/xue/software/trinityrnaseq-Trinity-v2.5.1//Analysis/DifferentialExpression/run_DE_analysis.pl --matrix ../liver.counts.matrix --method edgeR --samples_file ../liver_mvsf_samplefile_EdgeR.tsv
```


Output DE files were stored in: ` /home/xue/borealis_DE/edgeR_out`.

If I was to run kallisto without the Trinity scripts, I would do the following. But it is not neccessary to run it as it is essentially the same commands that Trinity was running. 
```
#building index for the assembled transcriptome
kallisto index -i /home/xue/borealis_DE/kallisto_index/Borealis_assembled_transcriptome_Trinity_kallisto_idx /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta

#quantification step of kallisto for each sample, can do quantification for one sample at a time
kallisto quant -i /home/xue/borealis_DE/kallisto_index/Borealis_assembled_transcriptome_Trinity_kallisto_idx -o /home/xue/borealis_DE/kallisto_index --plaintext -o female_rep7 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4082_girl_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4082_girl_liver_R2_scythe.fastq.gz)

#EdgeR

```


# Filter DE transcripts (FDR<0.05)
I filtered the edgeR output with a FDR treshold of 0.05 and select those transcripts that have at least fold change of 2 (or log2foldchange of 1). 
```
awk '($2 < -1||$2 >1) && $5<0.05  {print }' liver.counts.matrix.female_liver_vs_male_liver.edgeR.DE_results > liver_fdr005.tsv
```
A script were used to extract the sequence of significant DE transcripts and outputed the sequence to a fasta (.fa) file.
```
. get_trans_fdr005.sh
```
# Mapping DE X.borealis sequence against annotated X.laevis genome
### Mapping with Blastn
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

## Mapping with Star

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
Output file from STAR only contain mapping information for 8-30 mapping. I checked online and most people used STAR for mapping RNA-seq reads to genome but not mapping transcripts to genome. It might not be a tool that is capable for mapping transcripts to genome. Ben said he got ~80% mapping with STAR. I check with him and it seems like he had done read mapping but not transcripts mapping before. STAR doesn't seem like working for me, so switched to try other splice aware aligner such as GMAP or BLAT. 

### Mapping with GMAP
GMAP is another splice aware aligner and its user manual is here (http://research-pub.gene.com/gmap/src/README). It is already installed on the cluster (version 2012-04-27). 

```
#building local database index for *X. laevis* 9.2

#running gmap to map borealis transcriptome to laevis genome
gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -A -Z -f sampe /home/xue/laevis_dosageCompensated_gene/cap5n.fasta > cap5n_gmap.sam
#running gmap to map DE transcript to laevis genome
gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -A -Z -f sampe /home/xue/borealis_DE/liver_mvsf/filtered_edgeRout/liver_trans_fdr005_header.fa > /home/xue/borealis_DE/liver_mvsf/mapping_GMAP/liver_DE_gmap_out.sam
#running gmap to map laevis transcriptome to laevis genome


```

# Binning trancripts 
The Perl script that do the binning would: 
  - group overlapping transcripts into bins.  
  - count the number of transcripts in each bin. If 100k transcripts were mapped to chr8L, how many of them overlapped?
  - sum the expression level of transcripts in the same bin -> see if bins on chr8L non-recombinating region have different expression level than other bins (on chr8L non-sex-linked region and other autosomes)
  


# Compare to *X.laevis* chr8L
### *De novo* assembled *X.laevis* transcriptome
The path to *X.laevis* RNA-seq trimmed data and assembled transcriptome: /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Session-SRA. The path to the assembled transcriptome: /home/xue/borealis_DE/Session_Xl_DeNovo_Transcriptome  
```
/home/xue/software/trinityrnaseq-Trinity-v2.5.1/Trinity --seqType fq --left /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Session-SRA/*R1_scythe*  --right /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Session-SRA/*R2_scythe*  --CPU 20 --full_cleanup --max_memory 200G  --output /home/xue/borealis_DE/Session_Xl_DeNovo_Transcriptome/Xl_Transcript_TrinityOut
```

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
Error message:
```
Error: Encountered internal Bowtie 2 exception (#1)
Command: bowtie2-build --wrapper basic-0 --threads 20 -o 3 /home/xue/xborealis_alltogether/xb_liver_trinity/chrysalis/inchworm.K25.L25.DS.fa.min100 /home/xue/xborealis_alltogether/xb_liver_trinity/chrysalis/inchworm.K25.L25.DS.fa.min100
Error, cmd: bowtie2-build --threads 20 -o 3 /home/xue/xborealis_alltogether/xb_liver_trinity/chrysalis/inchworm.K25.L25.DS.fa.min100 /home/xue/xborealis_alltogether/xb_liver_trinity/chrysalis/inchworm.K25.L25.DS.fa.min100 1>/dev/null 2>tmp.39623.stderr died with ret 256 at /home/xue/software/trinityrnaseq-Trinity-v2.5.1/PerlLib/Pipeliner.pm line 166.
        Pipeliner::run(Pipeliner=HASH(0x1f440d0)) called at /home/xue/software/trinityrnaseq-Trinity-v2.5.1/Trinity line 1827
        main::run_chrysalis("/home/xue/xborealis_alltogether/xb_liver_trinity/inchworm.K25"..., "/home/xue/xborealis_alltogether/xb_liver_trinity/both.fa", 200, 500, undef, "/home/xue/xborealis_alltogether/xb_liver_trinity/both.fa", "/home/xue/xborealis_alltogether/xb_liver_trinity/both.fa") called at /home/xue/software/trinityrnaseq-Trinity-v2.5.1/Trinity line 1667
        main::run_Trinity() called at /home/xue/software/trinityrnaseq-Trinity-v2.5.1/Trinity line 1320
        eval {...} called at /home/xue/software/trinityrnaseq-Trinity-v2.5.1/Trinity line 1319

Trinity run failed. Must investigate error above.
```

# count the number of zero expresssion combinations
It was done with a perl script 'count_ZeroExpComb.pl'.
- all

|                |Female-0       |Female-1       |Female-2       |Female-3       |Female-4|
|---|---|---|---|---|---|
|Male-0        |5.93%  |3.22%  |1.09%  |0.46%  |0.14%|
|Male-1 |2.05%  |2.56%  |2.36%  |1.78%  |0.87%|
|Male-2 |0.68%  |2.06%  |3.43%  |4.27%  |3.44%|
|Male-3 |0.27%  |1.36%  |3.74%  |7.73%  |10.27%|
|Male-4 |0.09%  |0.57%  |2.65%  |8.96%  |30.03%|



- DE: total = 432 

|                |Female-0       |Female-1       |Female-2       |Female-3       |Female-4|
|---|---|---|---|---|---|
|Male-0 |3.92%  |0.46%  |1.15%  |1.61%  |14.98%|
|Male-1 |3.69%  |0      |0      |0.69%  |5.76%|
|Male-2 |3.00%  |0      |0      |0      |0|
|Male-3 |8.06%  |0.46%  |0      |0      |0|
|Male-4 |50.46% |5.76%  |0      |0      |0|


- on chr8L: total DE = 160


|		|Female-0	|Female-1	|Female-2	|Female-3	|Female-4|
|---|---|---|---|---|---|
|Male-0	|5.62%	|0	|0.62%	|1.88%	|7.50%|
|Male-1	|4.38%	|0	|0	|0.62%	|0.62%|
|Male-2	|5.00%	|0	|0	|0	|0|
|Male-3	|8.12%	|0	|0	|0	|0|
|Male-4	|64.38%	|1.25%	|0	|0	|0|

- on chr8L sex link region: total DE = 137


|		|Female-0	|Female-1	|Female-2	|Female-3	|Female-4|
|---|---|---|---|---|---|
|Male-0	|5.84%	|0	|0.73%	|1.46%	|8.03%|
|Male-1	|3.65%	|0	|0	|0.73%	|0.73%|
|Male-2	|5.11%	|0	|0	|0	|0|
|Male-3	|7.30%	|0	|0	|0	|0|
|Male-4	|64.96%	|1.46%	|0	|0	|0|
 
- on chr8L sex link region: total DE = 23

|		|Female-0	|Female-1	|Female-2	|Female-3	|Female-4|
|---|---|---|---|---|---|
|Male-0	|4.35%	|0	|0	|4.35%	|4.35%|
|Male-1	|8.70%	|0	|0	|0	|0|
|Male-2	|4.35%	|0	|0	|0	|0|
|Male-3	|13.04%	|0	|0	|0	|0|
|Male-4	|60.87%	|0	|0	|0	|0|

# things that we could do but not urgent
- repeat DE analysis with different softwware package (RESM, DESeq2 etc)
- repeat the pipeline with DE comparison pair (all male smaples vs all female samples, male gonad vs female gonad, liver vs gonad)

# Path to laevis RNAseq read data
```
/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Session-SRA
```

# Directionarity of expression evolution
## BM: Big Frame work
We want to examine the directionality of expression evolution of *X. borealis* SA gene by fitting them to BM model.

Phylogenetic tree:
- XB female/male expression ratio
- XL ortholog female/male expression ratio
- XT or XMT ortholog female/male expression ratio

We are going to test for 4 models:
-  all genes have 1 same sigma
- 1 rate for 8L sex linked region and 1 rate for the rest
- each gene have 2 sigma (one for XL and one for XB)
- all gene have 1 sigma but have different rate for different gene

We decided that PGLS might be a better model for this analysis. We are putting the BM model on hold.

# PGLS











