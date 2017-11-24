This is a summary of our conversations regarding assessing dosage compensation in *X.borealis* 

# Objective
- To assess dosage compensation in *X.borealis* with difference approaches.
- Quantify gene expression levels in tetraploid *X.b* in the absense of properly assembled genome

## *De novo* approach: mapping to assembled transcriptome
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

## Genome mapping approach 1: mapping to reference genome; allele specific expression
#### Method
- this is the approach that BenF is investigating
- build a reference *X.l* genome by pulling annotated genes from the *X.l* genome
- map assembled individual transcriptome to the genome with splicing aware aligner
- look for sex-linked SNPs that are heterozygous in female and homozygous in male
- assess the coverage for each SNP 
#### Advantage
- mapping to a reference genome would be a better alternative for mapping to a reference transcriptome.  
#### Caveats
- relies on the quality of *X.l* genome
- suggestion: get *X.l* coordinates from every *X.b* transcript and then collapse overlapping coordinates into a “gene" and sum the expression counts
  - problem: genes on the W chromosome is expected to have higher accumulation of mutation due to lack of recombination This will decrease mapping accuracy and underestimate gene expression level of genes on W chromosome compare to genes on Z chromosome or autosome. Thus, if expression level of W genes, we can't tell whether it is due to mapping error or true biological expression level. 
  - solution: find out what others are doing with estimating expression level of Z/W in the absence of reference genome

## Genome mapping approach 2: mapping to Unigene duplicates; locus expression
#### Questions
- are there different expression levels of genes (in the sex-link region) in males and in females? 
- are there more sex-biased genes in the sex-linked region than autosome regions?
- are there more sex-biased genes in the sex-linked region than the same region in the homeology chromosome?
- does the gene expression level of the sex-linked region (8L) in *X.b* higher compare to 8L in *X.l*?  
#### Method
- download genome of *X.laevis* duplicates from unigene database
- map the reads to the duplicates and estimate gene expression levels (read counts)
- DE analysis with RESEM (or others)
#### Advantages
- allow distinction between duplicates so that expression level of duplicates would not get overestimated 
#### Caveats
- singleton genes will not be included
- *X.l* duplicates might diverge enough that the *X.b* duplicates will only be mapped to one of the *X.l* duplicates
- the last update of the duplicates in the unigen database is 


