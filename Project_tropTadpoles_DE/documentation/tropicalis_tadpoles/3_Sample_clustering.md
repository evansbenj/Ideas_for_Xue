# Samples clustering - separate male and female samples
## Principle Component Analysis (PCA)
### Learning
I know nothing about PCA so I started this analysis by learning what it is and how to do it in R. 

#### What is PCA?

Reference
- StatQuest explained it very well, so I recommended watching this first:https://www.youtube.com/watch?v=FgakZw6K1QQ  
- Bio-Star Q&A: https://www.biostars.org/p/171925/
- nature methods: https://www.nature.com/articles/nmeth.4346
- Stackexchange Q&A: https://stats.stackexchange.com/questions/2691/making-sense-of-principal-component-analysis-eigenvectors-eigenvalues
- https://www.biostars.org/p/280615/

#### How to do PCA?
In the paper by Fehrmann et al (2015), their "preprocessing and aggregation of raw data were performed according to the robust multi-array average (RMA) algorithm in combination with quantile normalization".

Reference
- StatQuest explain the basic very well, I recommended watching this first: https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
- https://www.nature.com/articles/ng.3173#methods
- https://www.hindawi.com/journals/bmri/2018/2906292/
- https://www.datacamp.com/community/tutorials/pca-analysis-r
- https://www.biostars.org/p/282685/
- http://strata.uga.edu/8370/handouts/pcaTutorial.pdf
- http://uc-r.github.io/pca

#### PCA analysis
The script that I used to perform PCA is in the script folder with name `PCA_analysis_RNAseq_data.R`.

PCA analysis didn't really work in seperating our samples into two clusters. We are going to try more unsupervised clustering methods. A biostars Q&A (https://www.biostars.org/p/56897/) provided good resource to check out. 
## Kmeans
A phd student from Dr. Brian Golding's lab suggested kmean and set k=2 if we are looking for two group of clustering (males and females). 
#### What is Kmeans and how to do Kmeans?
- https://uc-r.github.io/kmeans_clustering
- https://2-bitbio.com/2017/10/clustering-rnaseq-data-using-k-means.html

## Identify cluster by examining expressions pattern of sex-related gene
### dmrt1
By theory, dmrt1 expression level should be different in tropicalis males and females. So Im going to blast dmrt1 sequence to the tropicalis transcriptome and found out which transcripts mapped to dmrt1 gene. Then I will look at the expression level of those transcripts. Our prediction is that they would have different expression level between sexes, which will allow us to identify the sex of the samples. We will then compare the list of males and females identified by dmrt1 epxression level to the list identified by clustering. The hope is that they should match to each other. If not, we need discuss about this further.
#### obtain sequence of dmrt1 genes in NCBI
I grabbed the mRNA and gene sequence of tropicalis dmrt1 gene from http://www.xenbase.org/gene/showgene.do?method=display&geneId=XB-GENE-852569. 
```
#path to the fasta file
/home/xue/tropicalis_gonad_transcriptome_Dec2018/analysis/samples_clustering_analysis/identify_sex_by_dmrt1/tropicalis_dmrt1_seq.fasta
```bash

**blastn dmrt1 sequence to the tropicalis transcriptome**
```
#blastn index tropicalis de novo transcriptome 
makeblastdb -in /home/xue/tropicalis_gonad_transcriptome_Dec2018/data/tropicali_gonad_transcriptome_trinityOut/tropicalis_transcriptome_build_dec2018/tropicalis_transcriptome_trinityOut.Trinity.fasta -dbtype nucl -out /home/xue/tropicalis_gonad_transcriptome_Dec2018/data/tropicali_gonad_transcriptome_trinityOut/tropicalis_transcriptome_build_dec2018/db_tropicalis_gonad_transcriptome_blastn/db_tropicalis_gonad_transcriptome

# blastn dmrt1 genes into tropicalis transcriptome
blastn -task blastn -db /home/xue/tropicalis_gonad_transcriptome_Dec2018/data/tropicali_gonad_transcriptome_trinityOut/tropicalis_transcriptome_build_dec2018/db_tropicalis_gonad_transcriptome_blastn/db_tropicalis_gonad_transcriptome -outfmt 6 -evalue 0.00005 -query  /home/xue/tropicalis_gonad_transcriptome_Dec2018/analysis/samples_clustering_analysis/identify_sex_by_dmrt1/tropicalis_dmrt1_seq.fasta -out /home/xue/tropicalis_gonad_transcriptome_Dec2018/analysis/samples_clustering_analysis/identify_sex_by_dmrt1/dmrt1_tropTrna_blastnOut.tsv
```

Then I ran the script `clustering_with_orthoGeneOnly.R` with blastn output `dmrt1_tropTrna_blastnOut.tsv`. Nothing really clutered.

### DE genes in borealis gonad (stage50)
Ben sent me the fasta file that contained sequence of DE gene (4 fold higher in female) in borealis stage 50 gonad. I uploaded it onto info and blastn it against the tropicalis transcriptome. 
```bash
blastn -task blastn -db /home/xue/tropicalis_gonad_transcriptome_Dec2018/data/tropicali_gonad_transcriptome_trinityOut/tropicalis_transcriptome_build_dec2018/db_tropicalis_gonad_transcriptome_blastn/db_tropicalis_gonad_transcriptome -outfmt 6 -evalue 0.00005 -query  /home/xue/tropicalis_gonad_transcriptome_Dec2018/analysis/samples_clustering_analysis/identify_sex_by_borealisGonad_DEgene/Femaletads_4X_maletads.fasta -out /home/xue/tropicalis_gonad_transcriptome_Dec2018/analysis/samples_clustering_analysis/identify_sex_by_borealisGonad_DEgene/borealisGonadDEgene_tropTrna_blastnOut.tsv
```
Then I ran the script `clustering_with_orthoGeneOnly.R` with blastn output `borealisGonadDEgene_tropTrna_blastnOut.tsv`. The R script would extract the transcript that mapped to those gene and perform PCA/Kmeans clustering with only expression value of those mapped transcripts.

The result from the above analysis produced the same result as the PCA and Kmeans analysis using the entire trop transcriptome expression data. Hence, Ben sent me the another fasta file that contain sequence of DE gene (rawcount cutoff of 4, and 4 fold higher in male or 4 folder higher in female) in borealis gonad at stage 50. I uploaded it onto info and blastn it against the tropicalis transcriptome.
```bash
#blastn; time cost: 44min
time blastn -task blastn -db /home/xue/tropicalis_gonad_transcriptome_Dec2018/data/tropicali_gonad_transcriptome_trinityOut/tropicalis_transcriptome_build_dec2018/db_tropicalis_gonad_transcriptome_blastn/db_tropicalis_gonad_transcriptome -outfmt 6 -evalue 0.00005 -query  /home/xue/tropicalis_gonad_transcriptome_Dec2018/analysis/samples_clustering_analysis/identify_sex_by_borealisGonad_DEgene/Sex_bias_4X_and_expression_4X.fasta -out /home/xue/tropicalis_gonad_transcriptome_Dec2018/analysis/samples_clustering_analysis/identify_sex_by_borealisGonad_DEgene/borealisGonadDEgene_4xBothSex_tropTran_blastnOut.tsv
```
Then I ran the script `clustering_with_orthoGeneOnly.R` with blastn output `borealisGonadDEgene_4xBothSex_tropTran_blastnOut.tsv`. PCA result remain the same - 3 different group of clutering, PC1 account for 12% of the variants (a little bit higher than the one above). Kmean produce a slightly different result. - with k=2, T11, T19, T20 clustered into a group; the rest clutered into a group 
- with k=3, the same 3 distinct group as before

## Using SNP
We have RADseq data for the east population that BenF did genotyping and found SNP that have a sex-specific pattern. We can mapped the RNAseq reads (East population) to the transcriptome and use the SNP information to identify the sex of the samples. 


