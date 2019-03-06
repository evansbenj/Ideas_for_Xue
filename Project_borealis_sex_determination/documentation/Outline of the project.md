# Outline

Based on conversation with Ben on Jan 8 2018, below is the overview of the bioinformatic steps involved in this project. 
- Build transcriptome
- Quality check using the genome and gtf file
- Collapse raw count to gene level with GTF 
  - we talked about collapsing the raw count to mRNA level, however, since some gene might have more splice variants than others, collapsing to mRNA levels could be cause uneven dilution of read counts.  
- Normalize expression value - normalize the fragment size: TMM 
- Normalize expression value - normalize with mean and standard deviation 
- TMM -> differential expression
  - DE for each developmental stage between male and female
- TMM + normalization with mean and standard deviation -> PCA
  - each dimension is one gene, and we are grouping them based on their expression value
  - one PCA for each developmental stage

