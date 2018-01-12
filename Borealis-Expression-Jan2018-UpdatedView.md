# *Xenopus Borealis* Gene Expression Study

There are (at least) 4 broad objectives of this project, outlined below. We have RNASeq data from parents and offspring for liver, and non-homologous sex specific tissues.

The main objective of this project is to test if newly evolved sex chromosomes have accumulated gene expression differences.

Some language to keep straight:

* With the transcriptome, we are really only looking at differential transcript expression, and not differential gene expression.


## 1. The Level of Sexual Antagonism (SA) in the Genome

**Objectives:**
	- the number of sex linked genes and their distribution in the genome
	- expression level of SA genes (i.e., are there huge expression differences between sexes?)
	- direction of expression differences (female up/down or male up/down), will need to leverage the *X. laevis* data to do this properly.

**Analysis Methods/Directions:**
	- Kallisto, full *de novo* approach
	- new Tuxedo package, full genomic approach

**Figures/Outputs/Analyses**
	- Volcano plot, colored by autosomes, 8L, and 8S
	- 2-fold expression, FDR.

**Things To Watch Out/Check For**
	-

## 2. Are sex chromosomes a hot-spot for SA? (related to \#1)

**Objectives**
	- percentage compared to the autosomes and 8S. Test if 8L has a disproportionate amount of SA genes
	- variation of expression level along 8L. Test if there is a peak somewhere (indicate earliest sex linked region?)
	- Can we tease apart alleles of the differentially expressed 8L genes. Is it an expression of the W copy in females high expression, both Z and W?
	- look at 8L in *X. laevis*, does it have more SA genes than other autosomes? (a speculation on my part that I am curious about.)

**Analysis Methods/Directions**
	- from \#1 above.

**Figures/Outputs/Analyses**
	- scatterplot 8L on x axis, expression level on y

**Things To Watch Out/Check For**
		- look for sex linked SNPs in 8S differentially expressed transcripts, could indicate mismapping.

## 3. Degeneration of sex chromosomes/allele specific expression (dosage compensation??)

**Objectives**
	- Are W and Z alleles expressed at different levels? All genes on chromosome, also look at differentially expressed genes specifically (see above).
	- Does the sex linked region have less expression overall in females? (Due to degeneration of W copy expression). If equal to males, could be dosage comp., but much more likely that W is just not degenerated. If we can show that females are expressing both genes then we know it's not upregulation of a Z copy.
	- Is 8L just suppressed expression in both sexes (balancing expression)? Could leverage laevis data to ask if females and males both have lower expression levels. Could compare to other the autosomes, ask if total median expression is lower.
	- Are there expression differences genome-wide, i.e., parent of origin effects on gene expression? (Het in one parent, hom in other, can track individual alleles).

**Analysis Methods/Directions**
	- mapping RNASeq to laevis genome
	- mapping to borealis superTranscriptome
	- currently just allele depths (normalization?)

**Figures/Outputs/Analyses**
	- box plot of expression levels of each allele and ZZ
	- probably linear mixed model (rand eff of pos so Z and W of each SNP are compared, could also include individual control for seq depth differences). Maybe a more simple paired-non parametic test (could even use median across individuals, but seq depth differences can affect this greatly).

**Things To Watch Out/Check For**
	- how do we know it's not just the W allele mapping more poorly? (limit analyses to alleles where W matches the ref -- took quick look, not looking favorable, think about it more).
	- look at expression levels of derived/ref alleles for non-8L sites. If mapping is problem, should see lower expression of derived allele.
	- look at boys ZZ sites, maternal SNPs can be tracked, expect not to see expression differences (unless some sort of parent of origin effect, like above). Look at effect of ZZ site being ref/diverged alleles.

## 4. What form the SA genes have? What does SA manifest as?

**Objectives**
	- trying to assess the mechanism of sexual antagonism, as opposed to just a description of its presence (which is sections \#1 & \#2)
	- splice variant differences between females and males? (i.e., both express the gene, just different forms)
	- whole gene expression differences (i.e., no expression of a gene in one sex)
	- sex specific subgenome expression, the roly of polyploidy in resolving sexual antagonism. Is it that females express an L copy and males express and S copy?

**Analysis Methods/Directions**
	- likely involve some genome mapping approaches
	- finding DE genes, and assessing if it is absent expression in one sex
	- location and overlap of DE genes (Xue's binning script will work for this)

**Figures/Outputs/Analyses**
	- figure or maybe table describing the characteristics of differentially expressed genes.
