# Role of splice variant in sex-biased transcript expression

## Objective
The purpose of this analysis is to see if splice variants contribute to the borealis DE distribution pattern.
There are two boroad area that we should look at. One is looking into the 432 DE borealis transcripts. The other one is look deeper into 160 collapsed gene. 

The things to check in the DE borealise transcripts that mapped to chr8L sex link region are:
- [ ] 1. how many transcript are lost when mapped to the annotated laevis gff genes?
- [ ] 2. how many of them are differential expressed when look at the transcript level but not DE when collasped into gene?
- [ ] 3. how many of them get merged when collasped into gene?

The things to check in the DE collapsed genes. 
- [ ] 1. how many of them contain mix of DE transcripts and non-DE transcripts? where are they?
- [ ] 2. how many of them contain non-DE transcripts only? - why do the gene become DE while it consist of non-DE transcript only?
- [ ] 3. 

The tings to check with non-DE collapsed genes.
- [ ] 1. how many of them contain DE transcripts but are non-DE collapsed genes? 
- [ ] 2. 


# What information do I have?

- list of DE transcripts
- list of DE collapsed genes
- list of collasped genes and the transcripts that mapped to each collapsed gene

# Look into DE borealis transcripts
### how many transcript are lost when mapped to the annotated laevis gff genes?
use R: read the bedfile of borealis transcript mapping to laevis gff genes; read the list of DE transcript that mapped to chr8L sex-link region. List/count the number of transcripts that are in both list. 

### how many of them are differential expressed when look at the transcript level but not DE when collasped into gene?
using R: 
