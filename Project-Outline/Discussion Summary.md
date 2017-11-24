This is a summary of our conversations regarding assessing dosage compensation in *X.borealis* 

# Objective
- To assess dosage compensation in *X.borealis* with difference approaches.
- Quantify gene expression levels in tetraploid *X.b* in the absense of properly assembled genome

# *De novo* approach: mapping to assembled transcriptome
#### Method
- mapping reads to assembled transcriptome that were assembled with reads from all *X.borealis* inididuals -> done
  - 4 comparisons: -> done
    - female tissues vs male tissues
    - female liver vs male liver
    - female gonad vs male gonad
    - liver vs gonad
- identify transcripts that have differentially expressed(DE) levels -> done
- blasting DE transcripts to *X.l* genome to identify genomic location 
#### Advantage
- will be able to capture genes that are present only in *X.b* but are missing in the *X.l* genome
#### Caveats
- large amount of transcripts in the assembled transcriptome could be allelic variants, splicing variants, and chimeric 
- we were talking about abandoming this approach; we think it is still valuable and maybe there are ways that we can improve it. I will look into literatures about estimating gene expression levels in polyploids (like plants) to see their way of doing it in the absence of properly assembled genome.  

# Genome mapping approach 
## Mapping to whole *X.l* genome
#### Method
- this is the approach that BenF is investigating
- map assembled individual transcriptome to the *x.l* genome with splicing aware aligner
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
- build a reference *X.l* genome by pulling annotated genes from the *X.l* genome using gff files
  - Alternative method for building a reference genome: get *X.l* coordinates from every *X.b* transcript and then collapse overlapping coordinates into a “gene" and sum the expression counts
- map the reads to the reference genome and estimate gene expression levels (read counts) with RESEM(or others)
- DE analysis with edgeR or DESeq2
#### Advantages
- allow distinction between duplicates so that expression level of duplicates would not get overestimated 


## Problem with genome mapping approach
- genes on the W chromosome is expected to have higher accumulation of mutation due to lack of recombination This will decrease mapping accuracy and underestimate gene expression level of genes on W chromosome compare to genes on Z chromosome or autosome. Thus, if expression level of W genes, we can't tell whether it is due to mapping error or true biological expression level. 
- *X.l* homeologs might diverge enough that both the *X.b* homeologs will be mapped to only one of the *X.l* hemeologs
- solution: find out what others are doing to estimate expression level of Z/W in the absence of reference genome


