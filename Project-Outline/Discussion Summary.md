This is a summary of our conversations regarding assessing dosage compensation in *X.borealis* 

# Objective
- To assess dosage compensation in *X.borealis* with difference approaches.
- Quantify gene expression levels in tetraploid *X.b* in the absense of high quality *X.b* genome

# *De novo* approach: mapping to assembled transcriptome
#### Method
- mapping reads to assembled transcriptome that was assembled with reads from all *X.borealis* inididuals -> done
  - 4 comparisons: -> done
    - female tissues vs male tissues
    - female liver vs male liver
    - female gonad vs male gonad
    - liver vs gonad
- identify transcripts that have differentially expressed(DE) levels: FC > 2, FDR <0.05;  -> done
- blast DE transcripts to *X.l* genome to identify chromosomal location 
#### Advantage
- would be able to capture genes that are present only in *X.b* or are missing in the *X.l* genome
#### Caveats
- large amount of transcripts in the assembled transcriptome could be allelic variants, splicing variants, and chimeric 
- we were talking about abandoning this approach; we think it is still valuable and maybe there are ways that we can improve it. I will look into literatures about estimating gene expression levels in polyploids (like plants) to see their way of doing it in the absence of properly assembled genome.  

# Genome mapping approach 
## Mapping to whole *X.l* genome
#### Method
- this is the approach that BenF is investigating
- map assembled individual transcriptome to the *X.l* genome with splicing aware aligner
- look for sex-linked SNPs that are heterozygous in female and homozygous in male
- assess the coverage for each SNP 

## Mapping to *X.l* reference genome
#### Questions
- are there different expression levels of genes (in the sex-link region) in males and in females? 
- are there more sex-biased genes in the sex-linked region than autosome regions?
- are there more sex-biased genes in the sex-linked region than the same region in the homeology chromosome?
- does the gene expression level of the sex-linked region (8L) in *X.b* higher compare to 8L in *X.l*?  
#### Method
- talked about download genome of *X.laevis* duplicates from unigene database -> abandoned because *X.l* files in Unigene database was previously updated in Mar 3 2013 
- build a *X.l* reference  genome by pulling annotated genes(exons) from the *X.l* genome based on gff3 files
  - Alternative method for building a reference genome: get *X.l* coordinates from every *X.b* transcript and then collapse overlapping coordinates into a “gene" and sum the expression counts
- map the assembled *X.b* transcriptome to the *X.l* reference  genome
- determine the consesus transcripts 
- Map read to the transcriptome and estimate gene expression levels with 
- DE analysis with edgeR or DESeq2

## Problem with genome mapping approach
- genes on the W chromosome are expected to accumulate more mutations due to lack of recombination. This will decrease mapping accuracy and underestimate gene expression level of genes on W chromosome compares to genes on Z chromosome or autosomes. Thus, if the analysis shows that W genes have lower expression level, we can't tell whether it is due to mapping errors or true biological variantions. 
- *X.l* homeologs might diverge enough that both *X.b* homeologs might be mapped to only one of the *X.l* hemeologs
- solution: find out what others are doing to estimate expression level of Z/W in the absence of reference genome


