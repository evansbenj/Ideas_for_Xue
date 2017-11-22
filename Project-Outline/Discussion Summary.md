This is a summary of our conversations regarding assessing dosage compensation in *X.borealis* 

# Objective
To assess dosage compensation in *X.borealis* through gene expression levels with difference approaches. 

# *De novo* approach: mapping to assembled transcriptome
### Method
- mapping reads to assembled transcriptome that were assembled with reads from all *X.borealis* inididuals
  - 4 comparisons:
    - female tissues vs male tissues
    - female liver vs male liver
    - female gonad vs male gonad
    - liver vs gonad
- identify transcripts that have differentially expressed(DE) levels
- blasting DE transcripts to *X.l* genome to identify genomic location 
### Advantage
- will be able to capture genes that are missing in the *X.l* genome
### Caveats
- large amount of transcripts in the assembled transcriptome could be allelic variants, splicing variants, and chimeric 
- we were talking about abandoming this approach; we think it is still valuable and maybe there are ways that we can improve it. I will look into literatures about DE in polyploids (like plants) to see their way of doing DE in the absence of properly assembled genome.  

# Genome mapping approach 1: mapping to reference genome; allele specific expression
### Method
- this is the approach that BenF is investigating
- build a reference *X.l* genome by pulling annotated genes from the *X.l* genome
- map assembled individual transcriptome to the genome with splicing aware aligner
- look for sex-linked SNPs that are heterozygous in female and homozygous in male
- assess the coverage for each SNP 
### Advantage
- 
### Caveats
- relies on the quality of *X.l* genome

# Genome mapping approach 2: mapping to Unigene duplicates; locus expression
### Method
- download genome of *X.laevis* duplicates from unigene database
- map the reads to the duplicates
- calculate the expression level of genes
### Advantages
- allow distinction between duplicates so that expression level of duplicates would not get overestimated 
### Caveats
- singleton genes will not be included
- *X.l* duplicates might diverge enough that the *X.b* duplicates will only be mapped to one of the *X.l* duplicates
