# Semi-automated differential expression workflow
Since I need to repeat the differential expression every time I make change to the up stream pipeline. I would like to make the DE process as automated as possible.

Steps:
1. read the input the file somewhat intelligently. or can have a user specific option
2. select the DE method: EdgeR or Deseq2
3. summarize the output of DE and output the summarize


# Reading input file
### example of input file
Example of a typical input file is a matrix with one culumn of gene/transcript ID and with columns of raw read count in each samples (ex, liver_female1, liver_female2). Below is an example of a matrix with gene ID and actual read count numbers.
```
Gene_id	female_rep1	female_rep2	female_rep3	female_rep4	male_rep1	male_rep2	male_rep3	male_rep4
MT:10711-11491	1394909	398621.9	240169.6	303280.5	377638	439351	310101.3	506429
MT:11562-11904	113992.269	34686.800000026	31881.4000006035	32997.1	51962.7661	187913.5732	111793.535	77881.472
MT:13855-15669	782.8771	264.4124	188.1878	63.48418	37449.08	525.469	206.5664	248.2572
MT:16249-17388	561066	236644	139493	151080	129554	333908	206843	252375
```
Below is an example of a matrix with transcript ID and actual read count numbers.
```
female_rep1	female_rep2	female_rep3	female_rep4	male_rep1	male_rep2	male_rep3	male_rep4
TRINITY_DN130196_c9_g1_i2	0	0	0	0	0	0	0	0
TRINITY_DN141981_c1_g4_i4	4.54574	0.000178133	9.02669	0	8.5702e-05	0	7.15789	9.39368
TRINITY_DN144064_c3_g1_i2	0	0	0	0	0	0	0	0
TRINITY_DN146823_c3_g1_i14	0	0	0	0	0	0	0	0
```
The reason that we should differential between matrix with gene ID and transcript ID is that with gene ID, we already know its genomic location, but with transcript ID, we need to read in the mapping file to find out the genomic locations.

### reading in the input file
Default is the matrix with transcript IDs. 
- [ ] input from user: path to the matrix
- [ ] getOption: matrix with transcript ID or gene ID? default is transcript ID
- [ ] a subroutine to read in mapping file **doesnt need to do at this step, can do it in the summaring step**

# Performing DE
### information needed
- read matrix: from above input
- filter: only sex-biased or inlcuded both sex-biased and sex-specific
- grouping: which columns belong to the same group; do it as what trinity do - a file that provide the grouping information
- method: edgeR or DESeq2? default is edgeR
- normalization
- output the result

# Summarizing DE result
### what type of summary do we want?
- summarize based on %DE per chr 
  -  output stat for each chr
  - output stat for categories; different categories for borealis, laevis, and tropicalis
    - borealis: chr8L sex link, chr8L non-sex link, chr8S, other chromosomes, scaffolds
    - laevis: chr2l sex link, chr2L non-sex link, chr2S, other chromosomes, scaffolds
    - tropicalis: chr7, other chromosomes, scaffolds
- 



