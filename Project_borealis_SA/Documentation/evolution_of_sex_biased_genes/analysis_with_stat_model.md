# Directionarity of expression evolution
## BM: Big Frame work
We want to examine the directionality of expression evolution of *X. borealis* SA gene by fitting them to BM model.

Phylogenetic tree:
- XB female/male expression ratio
- XL ortholog female/male expression ratio
- XT or XMT ortholog female/male expression ratio

We are going to test for 4 models:
-  all genes have 1 same sigma
- 1 rate for 8L sex linked region and 1 rate for the rest
- each gene have 2 sigma (one for XL and one for XB)
- all gene have 1 sigma but have different rate for different gene

We decided that PGLS might be a better model for this analysis. We are putting the BM model on hold.

# PGLS
### Input file
Input file format:
```
Xborealis_DE_Transcript_ID	borealis_F/M_ratio	laevis_F/M_ratio	tropicalis_F/M_ratio_randomNum
```


