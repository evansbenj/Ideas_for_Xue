# Differential exon usage

## DTU using Trinity wrapper
following tutorial: `https://github.com/trinityrnaseq/trinityrnaseq/wiki/DiffTranscriptUsage`

I tried to do it in info with the following:
```
/home/xue/software/trinityrnaseq-Trinity-v2.7.0-PRERELEASE/Analysis/SuperTranscripts/DTU/dexseq_wrapper.pl \
       --genes_fasta /home/xue/borealis_transcriptome/borealis_supertranscriptome_2018/supertranscriptome_dec2018/trinity_genes.fasta \
       --genes_gtf /home/xue/borealis_transcriptome/borealis_supertranscriptome_2018/supertranscriptome_dec2018/trinity_genes.gtf \
       --samples_file /home/xue/borealis_transcriptome/borealis_supertranscriptome_2018/borealis_superTrans_DTU/borealis_de_liver_samplefile.tsv \
       --out_prefix DTU --aligner STAR
```
It didn't work. 

I tried to run it on Graham. 
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3 
#SBATCH --time=01:12:00
#SBATCH --mem=100G
#SBATCH --job-name=borealis_supertranscriptome_dec2018
#SBATCH --account=def-ben

module load nixpkgs/16.09  gcc/7.3.0
module load nixpkgs/16.09  intel/2018.3
module load samtools/1.9
module load r/3.5.1
module load star/2.6.1c
module load python/3.7.0
module load trinity/2.8.4


/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/intel2016.4/trinity/2.8.4/trinityrnaseq-Trinity-v2.8.4/Analysis/SuperTranscripts/DTU/dexseq_wrapper.pl \
       --genes_fasta /home/songxy/projects/def-ben/songxy/borealis_transcriptome/supertranscriptome/trinity_genes.fasta \
       --genes_gtf /home/songxy/projects/def-ben/songxy/borealis_transcriptome/supertranscriptome/trinity_genes.gtf \
       --samples_file /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/borealis_liver_trimmed_reads_samples.txt \
       --CPU 20 \
       --out_prefix DTU --aligner STAR
```

```
#samples.txt

female female_rep1   /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE3897_mom_liver_R1_scythe.fastq.gz     /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE3897_mom_liver_R2_scythe.fastq.gz
female female_rep2   /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE4009_girl_liver_R1_scythe.fastq.gz     /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE4009_girl_liver_R2_scythe.fastq.gz
female female_rep3   /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE4072_girl_liver_R1_scythe.fastq.gz     /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE4072_girl_liver_R2_scythe.fastq.gz
female female_rep4   /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE4082_girl_liver_R1_scythe.fastq.gz     /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE4082_girl_liver_R2_scythe.fastq.gz
male male_rep1   /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE3896_dad_liver_R1_scythe.fastq.gz     /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE3896_dad_liver_R2_scythe.fastq.gz
male   male_rep2   /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE3929_boy_liver_R1_scythe.fastq.gz     /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE3929_boy_liver_R2_scythe.fastq.gz
male male_rep3   /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE4017_boy_liver_R1_scythe.fastq.gz     /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE4017_boy_liver_R2_scythe.fastq.gz
male male_rep4   /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE4039_boy_liver_R1_scythe.fastq.gz     /home/songxy/projects/def-ben/songxy/borealis_transcriptome/trimmed_reads/trimmed_reads_individual/BJE4039_boy_liver_R2_scythe.fastq.gz
```
It didn't run on Graham either.





# DTU with DEXseq
The usage of DEXSeq is following the DEXSeq vignett and also learned form the above link. 

Before running DEXSeq in R, count files needs to be reformatted to a style that DEXSeq can use. Since the transcripts of a gene might contain different exons, DEXSeq provide a python script that collapsed the exons information to define exon counting bins The first script converts the gtf file into a gff file with collapsed exon counting files, which should include the coordintate of exon bins. From DEXSeq vegnett section 2.4, this first script can be ran like this:

```bash
time python3 /home/xue/software/DEXSeq/inst/python_scripts/dexseq_prepare_annotation.py /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_dec2018/supertranscriptome/borealis_superTrans.gtf /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_dec2018/analysis/supertranscriptome/DTU/dexseq_count/borealis_SuperTrans_dexseq.gff
```

The second script generate actual counts. The input file includes:
- the gff file that was generated by the first script.  
- alignment file, which RNAseq read were mapped to super-transcriptome(or genome). The output bam/sam file need to be sorted by positions

The parameter that I used included:
- `-p yes` indicate that data is from a paired-end sequencing run
- `-s no` indicated that library is not strand specific
- `-f bam` indicated thtat the alignment file is in bam file format
- `-r pos` indicated that in the alignment file, paired-end data is sorted by
alignment position

```bash
time python3 /home/xue/software/DEXSeq/inst/python_scripts/dexseq_count.py -p yes -s no -f bam -r pos /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_dec2018/analysis/supertranscriptome/DTU/dexseq_count/borealis_SuperTrans_dexseq.gff /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_dec2018/analysis/supertranscriptome/mapping_reads_to_superTrans_star/BJE3896_dad_superTrans_star.bam /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_dec2018/analysis/supertranscriptome/DTU/dexseq_count/BJE3896_dad_count.txt
```
The second script `dexseq_count.py` ran forever - has been running for 2 months and still not done. This takes too long. Hence, I looked at alternative method to do this. 







## DTU using *featureCount* and DEXSeq - not done
This is done partially following the tutorial descripted in https://github.com/mikelove/rnaseqDTU/blob/master/vignettes/rnaseqDTU.Rmd.

The above script ran for two months and still going on. This is too slow. I came across this paper (https://www.biorxiv.org/content/biorxiv/early/2018/07/26/377762.full.pdf) which compare the speed, memory and accuracy about different DTU work flow. It concluded that dexseq_count.py has the slowest speed and featureCounts from the Rsubread package has the fastest speed. For the interest of time, I switch to the Rsubread-featureCount workflow even though it is not the standard workflow outlined in the DEXSeq manual. The featureCount workflow was used in this tutorial (https://github.com/mikelove/rnaseqDTU/blob/master/vignettes/rnaseqDTU.Rmd) 


Now, the count information is ready to be used by DESeq in a R script. Below is the R script. 

```{r}
library(rnaseqDTU)
csv.dir <- system.file("extdata", package="rnaseqDTU")

samps <- read.csv(file.path(csv.dir, "samples.csv"))
head(samps)
samps$condition <- factor(samps$condition)
table(samps$condition)
files <- file.path("/path/to/dir", samps$sample_id, "quant.sf")
names(files) <- samps$sample_id
head(files)

library(DEXSeq)
sample.data <- DRIMSeq::samples(d)
count.data <- round(as.matrix(counts(d)[,-c(1:2)]))
dxd <- DEXSeqDataSet(countData=count.data,
                     sampleData=sample.data,
                     design=~sample + exon + condition:exon,
                     featureID=counts(d)$feature_id,
                     groupID=counts(d)$gene_id)

system.time({
  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd, quiet=TRUE)
  dxd <- testForDEU(dxd, reducedModel=~sample + exon)
})

dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)
qval <- perGeneQValue(dxr)
dxr.g <- data.frame(gene=names(qval),qval)


columns <- c("featureID","groupID","pvalue")
dxr <- as.data.frame(dxr[,columns])
head(dxr)
```
### stageR 
```{r}
data(dex_tables)

library(stageR)
strp <- function(x) substr(x,1,15)
pConfirmation <- matrix(dxr$pvalue,ncol=1)
dimnames(pConfirmation) <- list(strp(dxr$featureID),"transcript")
pScreen <- qval
names(pScreen) <- strp(names(pScreen))
tx2gene <- as.data.frame(dxr[,c("featureID", "groupID")])
for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])

stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                      pScreenAdjusted=TRUE, tx2gene=tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05)
suppressWarnings({
  dex.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                 onlySignificantGenes=TRUE)
})
head(dex.padj)
```


How they do it in trinity
```
library(DEXSeq)
samples_info = read.table("$samples_table_file", header=T, row.names=1)
dxd = DEXSeqDataSetFromHTSeq(as.vector(samples_info\$counts_filename), sampleData=samples_info, design = ~ sample + exon + condition:exon, flattenedfile="$genes_gtf_file.dexseq.gff")
dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd )
plotDispEsts( dxd )
dxd = testForDEU( dxd )
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")
dxr1 = DEXSeqResults( dxd )
dxr1.sorted = dxr1[order(dxr1\$padj),]
save(list = ls(all=TRUE), file = "$out_prefix.Rdata")
write.table(dxr1.sorted, file="$out_prefix.results.dat", quote=F, sep="\t")

 pdf("$out_prefix.pdf")

 top_genes = unique(dxr1.sorted\$groupID[dxr1.sorted\$padj < 0.1 & ! is.na(dxr1.sorted\$padj)]);
 top_genes = top_genes[1:min($top_genes_plot, length(top_genes))];
 message("Top $top_genes_plot genes: (", paste(top_genes, collapse=','), ")") 
for (gene in top_genes) { 
   plotDEXSeq( dxr1 , gene, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 , expression=FALSE, norCounts=TRUE, splicing=TRUE, displayTranscripts=TRUE)
};

plotMA( dxr1, cex=0.8 )
```





## DTU using DRIMSeq - following rnaseqDTU tutorial 
https://github.com/mikelove/rnaseqDTU/blob/master/vignettes/rnaseqDTU.Rmd.


## Salmon
#### indexing
```
salmon index -t /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_dec2018/supertranscriptome/borealis_superTrans.fasta -i /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_dec2018/analysis/supertranscriptome/count_salmon/borealis_superTrans_index

```
#### quantifying  
bash script to run salmon for all the files
```bash
#!/bin/bash

index_dir=/home/xue/borealis_transcriptome/borealis_denovo_transcriptome_dec2018/analysis/supertranscriptome/count_salmon/borealis_superTrans_index

sample_dir=/home/xue/borealis_transcriptome/data/trimmed

out_dir=/home/xue/borealis_transcriptome/borealis_denovo_transcriptome_dec2018/analysis/supertranscriptome/count_salmon

cd /home/xue/borealis_transcriptome/data/trimmed

for i in *fastq.gz; do 
        sample_name=$(grep -o "BJE[0-9]*_[a-z]*" <(echo $i)) 
        
        salmon quant -i ${index_dir} -l A \
          -1 ${sample_dir}/${sample_name}_liver_R1_paired.fastq.gz \
          -2 ${sample_dir}/${sample_name}_liver_R2_paired.fastq.gz \
          -p 15 --validateMappings --rangeFactorizationBins 4 \
          --seqBias --gcBias \
          -o ${out_dir}/${sample_name}_quant
done
```
to run the bash script
```
~/script/run_salmon.sh
```




## DTU using DRIMSeq - following rnaseqDTU tutorial 
















