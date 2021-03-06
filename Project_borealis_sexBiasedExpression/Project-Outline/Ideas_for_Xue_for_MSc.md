# Ideas for Xue for M.Sc. work

Previous discussions with Xue about her upcomming M.Sc. work centered on using *X. laevis* to test for evidence of dosage compensation in amphibians. Ben Furman and I both realized that *X. borealis* would be a much better candidate to do this with because the region of suppressed recombination in *X. borealis* is much larger than in *X. laevis*.  So there should be many more genes to explore.

To this end we did RNAseq on liver tissue from a mother, a father, three daughters and three sons. These were multiplexed over at least two lanes (BenF can confirm). The point of doing this with a family is that we should be able to identify SNPs on the W and Z of the mother and distinguish the level of expression of each allele in the daughters and compare them individually to the Z alleles in the sons.

A major caveat of this project is that *X. boealis* may not have much divergence on it's sex chromosomes, which means *a priori* we do not expect there to be much expression differentiation.

Broad qustions I see (BenF):
1. Are the sex chromosomes a hot spot for sexually antagonism loci? (And, are the even any SA loci?)
2. Is there some degree of dosage compenstation? -- nobody has really assessed this in frogs, and again, *X. borealis* is not ideal, as there is minimal degeneration, so *a priori* there should be no reason for expression differences...but is still interesting to look at. 
3. Is the W allele more degraded than the Z allele -- i.e., more mutations?

# Collaborative nature of the project

This work will be a collaborative project that crucially involves BenF and both builds on and depends on the discoveries he made during his Ph.D. thesis work.  BenF will contribute to various aspects of this project as he wishes, and the division of labour will likely be a work in progress as the project goes on.

# Questions for Ben:
* Did we do the RNAseq on the smae individuals for which we have WGS data?
  - BenF says yes - cool!

* Did we do the RNAseq on only liver tissue or did we also include gonad?
  - BenF says no - cool!  We did liver on everyone and we did gonad on everyone except the mother (her RNAseq did not work)

* How many lanes did we do again?
  - BenF says 3 lane - each sample multiplexed across all three

* Did you back up the data on the lab tower?
  - BenF will do this.


# Main Directions

* We can test whether the overall expression level of the nonPAR alleles in the mother and daughters is similar to that in the father and sons.

* We can test in the mother, whether the expression level of the W allele differs from the Z allele, and whether each of these differs from the Z allele in the sons and father.

* We can test whether, in non-PAR sex linked genes, expression divergence between the sexes (the magnitude of expression divergence) tends to be more pronounced compared to the autosomes. The prediction if this region is a warehouse for genes with sexually antagonistic function is that there should be more sex biased expression divergence (even on Z alleles after controlling for copy number) than the autosomes.

* we can also do more detailed analysis of molecular evolution and splicing (which is the focus of the work the Xue did for her undergrad thesis) and divide this up in to genes with sex biased expression and non-sex biased expression.

# Pros and Cons

One of the neat aspects of doing this in *Xenopus borealis* is that all of these comparisons can be done in tandem with analysis of the homeologous genes on chromosome 8S, which is not sex-linked.

A challenge, I think, may be to ensure that we can distinguish between the homeologs on 8L and 8S with very high confidence.

Possible solution: when we need to distinguish, collect an *S. tropcialis* ortholog and the two *X. laevis* homeologs and build trees.


# Other

* Maybe get *X. laevis* transcriptomes to assess "proto-sex chromosome" dosage balances. I.e., what was the dosage level of the gene before being established as a sex chromosome? This gets at the idea of the "ancestral dosage" level.
	- that dosage review paper BenE sent said this "The hypothetically ideal experiment would be to compare the focal species with differentiated allosomes to a progenitor with undifferentiated proto-sex chromosomes."
	- with *X. laevis*, we could do exactly this. That paper suggests that given the wide variation in autosomal expression levels, comparing the sex chromosomes to the autosomes is not a good way to do the study. Instead, comparing it to an outgroup (like *X. laevis*) offers much more power and accuracy. For example, the sex chromosome may have higher expression than the autosomes, as when it was an autosome it had higher expression than the other autosomes. Therefore, if there was a dip in expression levels, we wouldn't detect it when comparing to the other autosomes but would detect it if comparing to the same chromosomes in an outgroup.

* We should spend some time digging around for available *X. laevis* transcriptomes (I need to do this for something else anyway). 

* Then again, we are also comparing males and females...so we'll be fine without *X. laevis* transcriptomes. They could just offer a different picture.


# Approaches for SA expression

I think it will be informative to take two approaches. One is a totally *de novo* differential expression analysis, where we find genes of interest, then map them back to the genome. The other approach would be to leverage the genome mapping from the start, and look for differences. These should converge on similar answers, of course, but may pick up some differences, i.e., some genes won't map to the genome (bad for method 2), or some transcripts may not be assembled correctly (bad for method 1). 
