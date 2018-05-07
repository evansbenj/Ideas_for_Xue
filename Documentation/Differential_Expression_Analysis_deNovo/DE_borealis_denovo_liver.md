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
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/util/abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix borealis_liver  --name_sample_by_basedir female_rep1/abundance.tsv female_rep2/abundance.tsv female_rep3/abundance.tsv female_rep4/abundance.tsv male_rep1/abundance.tsv male_rep2/abundance.tsv male_rep3/abundance.tsv male_rep4/abundance.tsv

time perl /home/xue/trinityrnaseq-2.2.0/util/abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix trans_counts  --name_sample_by_basedir BJE3896_dad_liver/abundance.tsv BJE3896_dad_testis/abundance.tsv BJE3897_mom_liver/abundance.tsv BJE3929_boy_liver/abundance.tsv BJE3929_boy_testis/abundance.tsv BJE4009_girl_liver/abundance.tsv BJE4009_girl_oviduct/abundance.tsv BJE4017_boy_liver/abundance.tsv BJE4017_boy_testis/abundance.tsv BJE4039_boy_liver/abundance.tsv BJE4039_boy_testis/abundance.tsv BJE4072_girl_liver/abundance.tsv BJE4072_girl_oviduct/abundance.tsv BJE4082_girl_liver/abundance.tsv BJE4082_girl_oviduct/abundance.tsv
```



Differential expression between two conditions (ex, male_liver samples vs female_liver samples) with EdgeR using Trinity scripts. Specifying EdgeR as the method as it is the recommended method :
```
time perl /home/xue/software/trinityrnaseq-2.2.0//Analysis/DifferentialExpression/run_DE_analysis.pl --matrix liver.counts.matrix --method edgeR --samples_file liver_mvsf_samplefile_EdgeR.tsv

time perl /home/xue/trinityrnaseq-2.2.0/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ../trans_counts.TMM.EXPR.matrix -P 1e-3 -C 2

#with a newer version of Trinity
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix /home/xue/borealis_DE/liver_mvsf/kallisto_out/borealis_liver.counts.matrix --method edgeR --samples_file /home/xue/borealis_DE/liver_mvsf/edgeR_out/borealis_liver_samples_files.tsv


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
#v4: total = 653
awk '($4 < -2||$4 >2) && $7<0.05  {print }' borealis_liver.counts.matrix.female_vs_male.edgeR.DE_results > liver_fdr005.tsv

#v2: 
awk '($2 < -1||$2 >1) && $5<0.05  {print }' liver.counts.matrix.female_liver_vs_male_liver.edgeR.DE_results > liver_fdr005.tsv
```
A script were used to extract the sequence of significant DE transcripts and outputed the sequence to a fasta (.fa) file.
```
. get_trans_fdr005.sh
```

# DE with Trinity_v4
```
#filter edgeR results and select transcripts that have |logfc|>1 and fdr<0.05; total number = 
awk '($4 < -2||$4 >2) && $7<0.05  {print }' borealis_liver.counts.matrix.female_vs_male.edgeR.DE_results > liver_edgeRout_de.tsv

#extract transcript sequences of de transcripts
perl ~/script/extract_sequence.pl liver_edgeRout_de.tsv /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta > liver_DE_seq.fa

#map those seq to laevis genome using blastn
blastn -task blastn -db ~/genome_data/laevis_genome/db_blastn_laevisGenome/Xl9_2_blastn_db -outfmt 6 -evalue 0.00005 -query  /home/xue/borealis_DE/liver_mvsf/edgeR_out/edgeR_trinity4/edgeR_liver_fvsm/liver_DE_seq.fa -out /home/xue/borealis_DE/liver_mvsf/edgeR_out/edgeR_trinity4/edgeR_liver_fvsm/liver_de_seq_blastout.tsv
```

# DE with Trinity_v51
```
#filter edgeR results and select transcripts that have |logfc|>1 and fdr<0.05; total number = 653
awk '($4 < -2||$4 >2) && $7<0.05  {print }' borealis_liver.counts.matrix.female_vs_male.edgeR.DE_results > liver_edgeRout_de.tsv

#extract transcript sequences of de transcripts
perl ~/script/extract_sequence.pl liver_edgeRout_de.tsv /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta > liver_DE_seq.fa
```


# Mapping DE X.borealis sequence against annotated X.laevis genome
There are a few options that we tried:
- blastn: not going to work well since it is not splice aware
- STAR: didn't work at first; BenF found that the package STAR-LONG works but dont have much documentation
- GMAP: tried with varies parameter; with DE transcripts, 93% mapping success; 
Others options:
- BLAT
- minimap2: very new (2018); from https://github.com/lh3/minimap2
- Exonerate: from https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate-manual


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
Output file from STAR only contain mapping information for 8-30 mapping. I checked online and most people used STAR for mapping RNA-seq reads to genome but not mapping transcripts to genome. It might not be a tool that is capable for mapping transcripts to genome. Ben said he got ~80% mapping with STAR. I check with him and it seems like he had done read mapping but not transcripts mapping before. STAR doesn't seem like working for me, so switched to try other splice aware aligner such as GMAP or BLAT. It seems like BenF get STAR to work with a package called STAR-Long

## Mapping with GMAP
GMAP is another splice aware aligner and its user manual is here (http://research-pub.gene.com/gmap/src/README). It is already installed on the cluster (version 2012-04-27). 

```
#building local database index for *X. laevis* 9.2

#running gmap to map borealis transcriptome to laevis genome (testing with subset of 100000 transcripts)
gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -A -B 5 -t 8 --cross-species -f samse /home/xue/borealis_DE/borealis_transcriptome/borealis_subset/borealis_subset_2.fasta | samtools view -S -b > /home/xue/borealis_DE/borealis_transcriptome/test/borealis_subset_2_gmap.bam

#running gmap to map DE transcript to laevis genome
gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -A -Z -f samse /home/xue/borealis_DE/liver_mvsf/filtered_edgeRout/liver_trans_fdr005_header.fa > /home/xue/borealis_DE/liver_mvsf/mapping_GMAP/liver_DE_gmap_out.sam

#testing to see which one will give me faster alignment, higher alignment success, and output in BAM file
gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -n 0 -A -B 5 -t 8 --cross-species -f samse /home/xue/borealis_DE/liver_mvsf/filtered_edgeRout/liver_trans_fdr005_header.fa | samtools view -S -b > sample.bam

gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -n 0 -A -B 5 -t 8 --cross-species -S /home/xue/borealis_DE/liver_mvsf/filtered_edgeRout/liver_trans_fdr005_header.fa > align_summary

gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -A -S -B 5 -t 8 --cross-species -S /home/xue/borealis_DE/liver_mvsf/filtered_edgeRout/liver_trans_fdr005_header.fa > alignment_only

#running gmap to map borealis transcriptome to laevis genome
gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -A -B 5 -t 8 --cross-species -f samse /home/xue/borealis_DE/borealis_transcriptome/borealis_subset/borealis_subset_17.fasta | samtools view -S -b > /home/xue/borealis_DE/borealis_transcriptome/test/borealis_subset_17_gmap.bam

```
When I tested with a subset of 100000 transcripts, it took 27719.79 seconds (7.7hrs, 3.61 queries/sec) and 43731/100000(43.7%) were unmapped. When I tested with DE transcripts, 31/432(7.2%) were not mapped. 
Below is a detailed note about parameter used in the above commands:
```
-D: path to the Genome directory
-d: name of the genome database
-A: Show alignments
-Z: Print output in compressed format -> probably dont want this if I want to pipe samtools to convert SAM to BAM
-B: batch mode; 5 is 
-f samse:specify output format; samse = SAM format (without setting paired_read bit); since the input is transcripts, they are single end. 
-t 6: number of thread that gmap can use; I indicated 6 thread/core
--cross-species: tell the program that it is a cross species alignment

```
Quick note about samtools: to count the number of unmapped transcripts
```
samtools view -f 0x4 sample.bam | wc -l
```

## Mapping with minimap

Install like this following instruction on https://github.com/lh3/minimap2
```
git clone https://github.com/lh3/minimap2
cd minimap2 && make
```
Indexing for laevis genome (/home/xue/genome_data/laevis_genome/db_minimap_laevis_92/XL9_2.mmi)
```
/home/xue/software/minimap2/minimap2 -d XL9_2.mmi ../XL9_2.fa
```
Trying mapping DE to laevis genome and see how long it will take
```
/home/xue/software/minimap2/minimap2 -cx asm20 --cs /home/xue/genome_data/laevis_genome/db_minimap_laevis_92/XL9_2.mmi /home/xue/borealis_DE/liver_mvsf/filtered_edgeRout/liver_trans_fdr005_header.fa > /home/xue/borealis_DE/liver_mvsf/mapping_minimap/liver_DE.paf
sort -k6,6 -k8,8n liver_DE.paf | /home/xue/software/minimap2/paftools.js call -f ecoli_ref.fa -L10000 -l1000 -> /home/xue/borealis_DE/liver_mvsf/mapping_minimap/liver_DE.vcf
```
Output of tester run: `/home/xue/borealis_DE/liver_mvsf/mapping_minimap/liver_DE.paf`

# Mapping borealis DE transcript to borealis genome
## mapping with Blastn to borealis genome_v1
Since we got the assembled borealis genome form Austin, I want to redo the mapping again
```
makeblastdb -in /home/xue/borealis_transcriptome/borealis_genome_v1/Xbo.v1.fa.gz -dbtype nucl -out /home/xue/borealis_transcriptome/borealis_genome_v1/db_borealis_genome_v1_blastn/borealis_genome_v1
```
## mapping with GMAP to borealis genome_v1
```
#building index
gmap_build -d borealis_genome_v1 -D /home/xue/borealis_transcriptome/borealis_genome_v1/db_borealis_genome_v1_gmap/ -g ../Xbo.v1.fa.gz 
#mapping
gmap gmap_build -d borealis_genome_v1 -D /home/xue/borealis_transcriptome/borealis_genome_v1/db_borealis_genome_v1_gmap/ -g ../Xbo.v1.fa.gz -A -Z -f samse /home/xue/borealis_DE/liver_mvsf/filtered_edgeRout/liver_trans_fdr005_header.fa > /home/xue/borealis_DE/de_sex_liver/mapping_GMAP_to_borealis_genome/liver_de_borealisGenome_gmap_out.sam
```


# count the number of zero expresssion combinations in DE transcripts
It was done with a perl script 'count_ZeroExpComb.pl'.
- all transcrips

|                |Female-0       |Female-1       |Female-2       |Female-3       |Female-4|
|---|---|---|---|---|---|
|Male-0        |5.93%  |3.22%  |1.09%  |0.46%  |0.14%|
|Male-1 |2.05%  |2.56%  |2.36%  |1.78%  |0.87%|
|Male-2 |0.68%  |2.06%  |3.43%  |4.27%  |3.44%|
|Male-3 |0.27%  |1.36%  |3.74%  |7.73%  |10.27%|
|Male-4 |0.09%  |0.57%  |2.65%  |8.96%  |30.03%|


- total DE: 432 

|                |Female-0       |Female-1       |Female-2       |Female-3       |Female-4|
|---|---|---|---|---|---|
|Male-0 |3.92%  |0.46%  |1.15%  |1.61%  |14.98%|
|Male-1 |3.69%  |0      |0      |0.69%  |5.76%|
|Male-2 |3.00%  |0      |0      |0      |0|
|Male-3 |8.06%  |0.46%  |0      |0      |0|
|Male-4 |50.46% |5.76%  |0      |0      |0|


- DE on chr8L: 160


|		|Female-0	|Female-1	|Female-2	|Female-3	|Female-4|
|---|---|---|---|---|---|
|Male-0	|5.62%	|0	|0.62%	|1.88%	|7.50%|
|Male-1	|4.38%	|0	|0	|0.62%	|0.62%|
|Male-2	|5.00%	|0	|0	|0	|0|
|Male-3	|8.12%	|0	|0	|0	|0|
|Male-4	|64.38%	|1.25%	|0	|0	|0|

- DE on chr8L sex link region: 137


|		|Female-0	|Female-1	|Female-2	|Female-3	|Female-4|
|---|---|---|---|---|---|
|Male-0	|5.84%	|0	|0.73%	|1.46%	|8.03%|
|Male-1	|3.65%	|0	|0	|0.73%	|0.73%|
|Male-2	|5.11%	|0	|0	|0	|0|
|Male-3	|7.30%	|0	|0	|0	|0|
|Male-4	|64.96%	|1.46%	|0	|0	|0|
 
- on chr8L non-sex link region: total DE = 23

|		|Female-0	|Female-1	|Female-2	|Female-3	|Female-4|
|---|---|---|---|---|---|
|Male-0	|4.35%	|0	|0	|4.35%	|4.35%|
|Male-1	|8.70%	|0	|0	|0	|0|
|Male-2	|4.35%	|0	|0	|0	|0|
|Male-3	|13.04%	|0	|0	|0	|0|
|Male-4	|60.87%	|0	|0	|0	|0|

