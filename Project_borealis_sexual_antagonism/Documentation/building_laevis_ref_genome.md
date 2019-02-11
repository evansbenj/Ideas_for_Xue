This page will document the commands and progress of qualifying gene expression level of *Xenopus borealis* through a reference genome mapping approach.

# Building *X.laevis* reference genome
gff3 file includes coordinate information of genes, RNAs, exons, and CDS. I will use X.laevis gff3 to extract transcript sequences from X.laevis genome. The newest version(9.2) of X.laevis gff3 and genome were downloaded form Xenbase:
```
wget ftp://ftp.xenbase.org/pub/Genomics/JGI/Xenla9.2/XENLA_9.2_Xenbase.gff3
wget ftp://ftp.xenbase.org/pub/Genomics/JGI/Xenla9.2/XL9_2.fa.gz
```
The format of gff3 file:
> Line 0: GFF version
>
> line 1: define the boundaries of region being annotated (gene name, start, end)
> 
> Line 2: gene annotations start, each row contains information for the the following columns (separated by tap):  
> 1. **seqid**:the name of the sequence(ex.chromosome/scafold) on which the feature exists
> 2. **source**: the algorithm or operating procedure that generate the featuer. Ex,"Genescan", "Genebank"
> 3. **type**: sequence ontology or an SO accession number; ex, "CDS", "stop codon", "exon"
> 4. **start**: the start coordinates of the sequence in the chromosome or scafold
> 5. **end**: the end coordinate of the sequence in the chromosome or scafold
> 6. **score**: the e-value
> 7. **strand**: +/- or . for unknown/not stranded
> 8. **phase**: reading frame indicated by 0-2 if the type is "CDS"
> 9. **attributes**: a list of information related to the seqid; separated by ";". Information such as ID, alias, parents, name, gap

### Extract transcripts 
X.laevis gff3 file includes splice vairiants of a gene as differernt mRNAs. Splice variants have different combination of exons from the same gene. There are two ways that I can extract the transcripts. One way is to extract all mRNAs and preverse the splicing variants of each gene. Another way is that I can extract exons from mRNAs (of the same gene) and collapse them into a transcript. I will lose splice variants information by doing latter but I will be able to capture splice variants that are present in X.borealis and absence in X.laevis. I talked to BenF about this. We didn't come to a conclusion on which way is better. I decide to proceed with the former for now (**will talke to BenE about this when skype**).       
I am trying to extract mRNA sequences based on annotation information in the X.laevis gff3 file. Here are some ways that others tried:
```
http://cole-trapnell-lab.github.io/cufflinks/file_formats/
https://www.biostars.org/p/196047/
https://deweylab.github.io/RSEM/rsem-prepare-reference.html
http://www.gettinggeneticsdone.com/2014/03/visualize-coverage-exome-targeted-ngs-bedtools.html
https://www.biostars.org/p/46281/
http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/extractfeat.html
```
**method 1: use gffread command**

I ran gffread as below. gffread extracted mRNA based on CDS coordinates and not exon coordiantes (about gffread: http://ccb.jhu.edu/software/stringtie/gff.shtml). Since CDS sequences are shorter than exons sequences and don't include sequences of 5'UTR and 3'UTR, the UTR sequences will be lost. This might not be ideal since RNA-seq does sequence UTR regions, which will be present in the transcriptome. Actually, gffread does extract just the exons by using flage -w (-c for CDS).  
```
gffread/gffread/gffread XENLA_9.2_Xenbase.gff3 -g XL9_2.fa -w xl_mRNA.fasta 
```
Just as a cautious test, I removed all the CDS information in the gff3 and then do the sequence extraction again. I compared the two files and it seems like they are the same.  
```
awk '$3 != "CDS"{print}' XENLA_9.2_Xenbase.gff3 > XENLA_noCDS.gff3
gffread/gffread/gffread XENLA_noCDS.gff3 -g XL9_2.fa -w xl_mRNA_test.fasta 
```

**method 2: Bedops + bedtools**

About Bedop: http://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/gff2bed.html

About Bentools: http://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html

I did some initial filtering and create a gff file that contains only the exons.
```
grep -w "exon" XENLA_9.2_Xenbase.gff3 > XENLA_exon.gff3
```
then pulled the sequence of exons using bedtools
```
bedtools getfasta -fi XL9_2.fa -bed XENLA_exon.gff3 -fo exons.fa
```
Acutally, gffread do extract transcripts with exon coordinates. gffread doesn't need additional filtering step thus, I am not going to use the bedops+bedtools method. However this method might be useful in extracting individual gene sequence with all exons (no splice variants). 

# Mapping transcriptome to referencce genome
Before starting the alignment, database pre/indexing was done with the following:
```
bwa index -p XLmrna_bwa_db -a is xl_mRNA.fasta
```
The alignment was done with BWA as below. It took 111 minutes to run. The default output of BWA is SAM file (SAM format detail: https://samtools.github.io/hts-specs/SAMv1.pdf).  Output files was kept in bam file (to save storage space).  
```
time bwa mem XLmrna_bwa_db /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta |samtools view -b > bwa_XBtoXL_output.bam 
```
In a SAM file, the second column of each row is the flag column, which indicates the status of the alignment (ex, mapped or unmapped). If the flag is 4, it means that the sequence is unmapped. I filtered out the transcripts that is unmapped.
```
samtools view -b -F 4 bwa_XBtoXL_output.bam |samtools sort -o XBtoXL_sorted.bam
```
Then I did some more filtering and genotyping. Filter out: multimapping transcripts, mapping quality. genotyping: BCF tool: skip indel/SNP sites (-V indels), alternative modelfor multiallelic and rare-variant calling(-m), output format is compressed VCF (-O z) 
```
samtools mpileup -ugf xl_mRNA.fasta -t DP,AD -s XBtoXL_sorted.bam |bcftools call -V indels --format-fields GQ -m -O z -o XBtoXL_consensus.vcf
```
The total number of transcript IDs is 1629500 before filtering and is 1062075 after filtering. 

As a note for future use______________________________________________________________

### Differential expression

discussion about FPKM vs counts vs RPKM
```
http://www.cureffi.org/2013/09/12/counts-vs-fpkms-in-rna-seq/
```

### other aligner
- Stampy: http://www.well.ox.ac.uk/~gerton/README.txt
- lastz: http://www.bx.psu.edu/~rsharris/lastz/README.lastz-1.04.00.html
- mummer: http://mummer.sourceforge.net/

