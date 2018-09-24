# Callapsing transcripts
Collapsing the 1.6 million transcripts to a number that is more biologically realistic (i.e. remove all of the false transcripts due to sequence error, polymorphism, etc..).

The callapsing of transcripts was done using a perl script. The script can be found in the script folder on info (path: /home/xue/script/identify_genomic_location_for_transcripts.pl) . Here I will illustrate the basic sctructure of the script. 

 - This perl script will take in the mapping file, in which the transcriptome was mapped to a genome and will provide the genomic location of mapping.
 - With the aid of the laevis gff file, which contain genomic coordinates of all annotated genes, the script can figure out which gene the transcripts mapped to based on genomic coordinate.
 - Afterward, the script will read in the expression value of each transcript (output of kallisto) and sum the expression value of transcripts that mapped to the same gene.
 - Lastly, it will output the callapsed expression value for each gene. An example of the output file is at the end of this document. 


#### generate gene only gff file 
extract gene information from gff file
```
awk '$3 == "gene"{print}' XENLA_9.2_Xenbase.gff3 > /home/xue/genome_data/laevis_genome/gff3_xl9_2/XENLA_gene.gff3
```

##### generate the alignment file
mapping transcriptome to genome
```
time gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -A -B 5 -t 8 -f samse --cross-species /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/subset_1.fasta | samtools view -S -b > /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/borealis_denovoT_laevisV92_genome_gmap_subset1.bam
```
filter the mapping read remove the unmapped reads, which is indicated by flag 0x04. Then I piled the filtered output to bedtools, which extracts alignment coordinates for transcripts based on their CIGAR strings;
```
samtools view -F 0x04 -b borealis_denovoT_laevisV92_genome_gmap.bam | bedtools bamtobed -i > borealis_denovoT_laevisV92_genome_gmap_bedfile.bed
```
filter the bedfile and keep only the hit with the highest mapping quality score; if there are two hit with the same mapping quality score, keep the first one only;
```
perl ~/script/orthologs_identification_filterBedfile.pl borealis_denovoT_laevisV92_genome_gmap_bedfile.bed > borealis_de_laevisV92_genome_gmap_bedfile_filtered.tsv
```
##### expression matrix 
```
#Indexing the transcriptomes:
kallisto index -i /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta.kallisto_idx /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta
# Transcript abundance quantification:
kallisto quant -i /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta.kallisto_idx  -o female_rep7 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4082_girl_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4082_girl_liver_R2_scythe.fastq.gz)
# compile the read count for each samples into a matrix; used a script included in the Trinity RNAseq analysis package
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/util/abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix borealis_liver  --name_sample_by_basedir female_rep1/abundance.tsv female_rep2/abundance.tsv female_rep3/abundance.tsv female_rep4/abundance.tsv male_rep1/abundance.tsv male_rep2/abundance.tsv male_rep3/abundance.tsv male_rep4/abundance.tsv
```

#### Path to all my input files
```
#laevis gff
#####/home/xue/genome_data/laevis_genome/gff3_xl9_2/XENLA_gene.gff3
#alignment files: borealis, laevis, tropicals
### 
/home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/borealis_denovoT_laevisV92_genome_gmap_filtered.bed
/home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_laevisGenomeApproach/laevis_denovoT_laevisV92_genome_gmap_bedfile_filtered.bed
/home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_laevisGenomeApproach/tropicalis_denovoT_laevisV92_genome_gmap_bedfile_filtered.bed

#expression matrix: borealis, laevis, tropicalis
/home/xue/borealis_DE/de_sex_liver/kallisto_out/borealis_liver.counts.matrix
/home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/kallisto_out/laevis_denovo.counts.matrix
/home/xue/tropicalis_transcriptome/tropicalis_denovo_transcriptome/tropicalis_kallisto_denovo/tropicalis_denovo.counts.matrix
```

Path to all the output files
```
/home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/sum_expression/borealis_expression_matrix_per_gene.tsv

/home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/sum_expression/laevis_expression_matrix_per_gene.tsv

/home/xue/tropicalis_transcriptome/tropicalis_denovo_transcriptome/sum_expression/tropicalis_expression_matrix_per_gene.tsv
```

My usage of the perl script
```
time perl ~/script/identify_genomic_location_for_transcripts.pl

time perl ~/script/identify_genomic_location_for_transcripts.pl --align /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/borealis_denovoT_laevisV92_genome_gmap_filtered.bed --matrix /home/xue/borealis_DE/de_sex_liver/kallisto_out/borealis_liver.counts.matrix > borealis_expression_matrix_per_gene_allowOverlap.tsv

time perl ~/script/identify_genomic_location_for_transcripts.pl --align /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_laevisGenomeApproach/laevis_denovoT_laevisV92_genome_gmap_bedfile_filtered.bed --matrix /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/kallisto_out/laevis_denovo.counts.matrix > /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/sum_expression/laevis_expression_matrix_per_gene_allowOverlap.tsv

time perl ~/script/identify_genomic_location_for_transcripts.pl --align /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_laevisGenomeApproach/tropicalis_denovoT_laevisV92_genome_gmap_bedfile_filtered.bed --matrix /home/xue/tropicalis_transcriptome/tropicalis_denovo_transcriptome/tropicalis_kallisto_denovo/tropicalis_denovo.counts.matrix > /home/xue/tropicalis_transcriptome/tropicalis_denovo_transcriptome/sum_expression/tropicalis_expression_matrix_per_gene_allowOverlap.tsv
```

An example of the output format 
```
        female_rep1     female_rep2     female_rep3     female_rep4     male_rep1     male_rep2     male_rep3     male_rep4
MT:10711-11491  973972  299235  179192  227215  276145  327805  229833  375904
MT:16249-17388  561066  236644  139493  151080  129554  333908  206843  252375
Scaffold1000:541-2366   0       0       1       0       0       0       0       0
Scaffold100:110725-181202       0       11      9       6       4       10      21      1
Scaffold100:186956-193508       7       1       6       2       2.41232 1       0       5.76162
Scaffold100:209407-303984       43      105.5498        58      73      98      33.20251        136     40.0171



```
### DE of collapsed transcripts


DE with edgeR and DESeq2
```
time perl /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix borealis_expression_per_gene.matrix --method edgeR --samples_file /home/xue/borealis_DE/de_sex_liver/edgeR_out/borealis_liver_samples_files.tsv

time perl /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix borealis_expression_per_gene.matrix --method DESeq2 --samples_file /home/xue/borealis_DE/de_sex_liver/edgeR_out/borealis_liver_samples_files.tsv

```
filter by FDR
```
awk ' $10<0.05  {print }' borealis_expression_per_gene.matrix.female_vs_male.DESeq2.DE_results > liver_deseq2_fdr005.tsv
awk '($6 < -2||$6 >2) && $10<0.05  {print }' borealis_expression_per_gene.matrix.female_vs_male.DESeq2.DE_results > liver_logfc2_fdr005.tsv

awk ' $7<0.45  {print }' borealis_expression_per_gene.matrix.female_vs_male.edgeR.DE_results > liver_edgeR_fdr045.tsv
awk '($4 < -2||$4 >2) && $7<0.05  {print }' borealis_expression_per_gene.matrix.female_vs_male.edgeR.DE_results > liver_edgeR_logfc2_fdr005.tsv


```
summarized the filtered result
```
perl ~/script/summarized_DE_result.pl liver_edgeR_fdr005.tsv

```


