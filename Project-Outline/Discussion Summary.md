This is a summary of our conversations regarding the *X.borealis* dosage compensation project

# Objective
To assess dosage compensation in *X.borealis* using gene expression level with difference approaches. 

# De novo approach: mapping to assembled transcriptome
### Method
- mapping reads to assembled transcriptome that were assembled with reads from all X.l inididuals
-  
### Advantage
- will be able to capture genes that is missing in the X.l genome
### Caveats
- large amount of transcripts in the assembled transcriptome could be allelic variants, splicing variants, and chimeric 

# Genome mapping approach 2: mapping to reference genome
### Method
- build a reference genome by pulling annotated genes from the genome
- 
### Advantage
- 
### Caveats
- 

# Genome mapping approach 1: mapping to Unigene duplicates
### Method
- download genome of *X.laevis* duplicates from unigene database
- map the reads to the duplicates
- calculate the expression level of genes
### Advantages
-
### Caveats
- singleton genes will not be included
- X.l duplicates might be diverge enough that the X.b duplicates will only be mapped to one of the X.l duplicates
- 

