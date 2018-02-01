# *Xenopus Borealis* Gene Expression Study

There are (at least) 4 broad objectives of this project, outlined below. We have RNASeq data from parents and offspring for adult liver, and non-homologous sex specific tissues. We soon will have comparable data from adult liver of 3 male and 3 female XL.

The main objective of this project is to test if newly evolved sex chromosomes have accumulated gene expression differences. This is a test that recombination suppression is favored by natural selection as a way to resolve sexual antagonism. Also want to look at allelic specific expression to test for degeneration of sex chromosomes (maybe a little dosage comp stuff, if we find differences).

Some language to keep straight:

* With the transcriptome, we are really only looking at differential transcript expression, and not differential gene expression.
* Once we have expression from the same gene in two species, and potentially two paralogs from one or both, we can make inferences about the directionality of expression evolution (e.g. this gene was upregulated or downregulated in females, etc).


## 1. The Level of Sexual Antagonism (SA) in the Genome

**Objectives:**
* the number of sex linked genes and their distribution in the genome

* expression level of SA genes (i.e., are there huge expression differences between sexes?)

* direction of expression differences (female up/down or male up/down), will need to leverage the *X. laevis* data to do this properly.

**Analysis Methods/Directions:**
* Kallisto, full *de novo* approach
** I'd like more details on how this approach works
** Also need to look into syncing data from XL and XB
*** We plan to generate expression levels from each species using transcriptome assemblies from each species, and then sync the noramlized data after this in a way that is not influenced by the sequencing depth of each dataset (possibly using reads per million).

* new Tuxedo package, full genomic approach

**Figures/Outputs/Analyses**
* Volcano plot, colored by autosomes, 8L, and 8S

* 2-fold expression, FDR.

**Things To Watch Out/Check For**

*

## 2. Are sex chromosomes a hot-spot for SA? (related to \#1)

**Objectives**

* percentage compared to the autosomes and 8S. Test if 8L has a disproportionate amount of SA genes

* variation of expression level, also variation in level of sex biased expression along 8L and compared to autosomes. Test if there is a peak somewhere (indicate earliest sex linked region?)

* Can we tease apart alleles of the differentially expressed 8L genes. Is it an expression of the W copy in females high expression, both Z and W?

* look at 8L in *X. laevis*, does it have more SA genes than other autosomes? (a speculation on my part that I am curious about.)  And this is predicted by thory - genomic regions with more sexual antagonism presumably should be favored by natural selection to become new sex chromosomes.

**Analysis Methods/Directions**

* from \#1 above.

**Figures/Outputs/Analyses**

* scatterplot 8L on x axis, expression level on y

**Things To Watch Out/Check For**

* look for sex linked SNPs in 8S differentially expressed transcripts, could indicate mismapping.

## 3. Degeneration of sex chromosomes/allele specific expression (dosage compensation??)

**Objectives**

* Are W and Z alleles expressed at different levels? All genes on chromosome, also look at differentially expressed genes specifically (see above).

* Does the sex linked region have less expression overall in females? (Due to degeneration of W copy expression). If equal to males, could be dosage comp., but much more likely that W is just not degenerated. If we can show that females are expressing both genes then we know it's not upregulation of a Z copy.

* Is 8L just suppressed expression in both sexes (balancing expression)? Could leverage laevis data to ask if females and males both have lower expression levels. Could compare to other the autosomes, ask if total median expression is lower.

* Are there expression differences genome-wide, i.e., parent of origin effects on gene expression? (Het in one parent, hom in other, can track individual alleles).

**Analysis Methods/Directions**

* mapping RNASeq to laevis genome

* mapping to borealis superTranscriptome

* currently just allele depths (normalization?)

**Figures/Outputs/Analyses**

* box plot of expression levels of each allele and ZZ

* probably linear mixed model (rand eff of pos so Z and W of each SNP are compared, could also include individual control for seq depth differences). Maybe a more simple paired-non parametic test (could even use median across individuals, but seq depth differences can affect this greatly).

**Things To Watch Out/Check For**

* how do we know it's not just the W allele mapping more poorly? (limit analyses to alleles where W matches the ref -- took quick look, not looking favorable, think about it more).

* look at expression levels of derived/ref alleles for non-8L sites. If mapping is problem, should see lower expression of derived allele.

* look at boys ZZ sites, maternal SNPs can be tracked, expect not to see expression differences (unless some sort of parent of origin effect, like above). Look at effect of ZZ site being ref/diverged alleles.

## 4. What form the SA genes have? What does SA manifest as?

**Objectives**

* trying to assess the mechanism of sexual antagonism, as opposed to just a description of its presence (which is sections \#1 & \#2)

* splice variant differences between females and males? (i.e., both express the gene, just different forms). Theory predicts less differential expression of splice variants in the sex-linked region because SA is already resolved by recombination suppression.

* whole gene expression differences (i.e., no expression of a gene in one sex). Expected to be more frequent in suppresed recombination region 

* sex specific subgenome expression, the roly of polyploidy in resolving sexual antagonism. Is it that females express an L copy and males express and S copy?  Interesting idea - unlikely I think but worth checking.

* If we find interesting genes in *X. borealis*, we could look at their expression across *X. laevis* tissues (Sessions et al. data) to see if they have broad or narrow expression profiles. 

* With interesting genes (IDed by liver transcriptomes), we can look at expression in sex specific tissues to see if they are high in one and absent in the other -- suggestive of sex specific function (we expect female specific genes to accumulate on the W).  Possible, but not so interesting?

**Analysis Methods/Directions**

* likely involve some genome mapping approaches

* finding DE genes, and assessing if it is absent expression in one sex

* location and overlap of DE genes (Xue's binning script will work for this)

**Figures/Outputs/Analyses**

* figure or maybe table describing the characteristics of differentially expressed genes.

## 5. Other sources of information

Later possibly compare rates of evolution (in nucleotides or dN/dS) of genes with and without sex biased expression patterns in and out of suppressed recombi region.



## Learning from the pro (beautiful!)

<img src="https://github.com/evansbenj/Ideas_for_Xue/blob/master/misc/benevans-style.jpg" width="600" />
