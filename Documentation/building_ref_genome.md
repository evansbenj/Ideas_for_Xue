This page will document the commands and progress of qualifying gene expression level of *Xenopus borealis* through a reference genome mapping approach.

# Building *X.laevis* reference genome
gff3 file of X.laevis was downloaded form Xenbase:
```
wget ftp://ftp.xenbase.org/pub/Genomics/JGI/Xenla9.2/XENLA_9.2_Xenbase.gff3
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

### Extract transcripts using bedtool
X.laevis gff3 file includes splice vairiants of a gene as differernt mRNAs. Splice variants have different combination of exons from the same gene. There are two ways that I can extract the transcripts. One way is to extract all mRNAs and preverse the splicing variants of each gene. Another way is that I can extract exons from mRNAs (of the same gene) and collapse them into a transcript. I will lose splice variants information by doing latter but I will be able to capture splice variants that are present in X.borealis and absence in X.laevis. I talked to BenF about this. We didn't come to a conclusion on which way is better. I decide to proceed with the former for now (**will talke to BenE about this when skype**).       
I am trying to extract mRNA sequences based on annotation information in the X.laevis gff3 file. Here are some ways that others tried:
```
http://cole-trapnell-lab.github.io/cufflinks/file_formats/
https://www.biostars.org/p/196047/
```
Way 1: use gffread command
I ran gffread as below. gffread extracted mRNA based on CDS coordinates and not exon coordiantes. Since CDS sequences are shorter than exons sequences and don't include sequences of 5'UTR and 3'UTR, the UTR sequences will be lost. This might not be ideasl since RNA-seq does sequence UTR regions, which will be present in the transcriptome.  


Way2: Bedpods + bedtools
I did some initial filtering and create a gff file that contains only the exons.
```
grep -w "exon" XENLA_9.2_Xenbase.gff3 > XENLA_exon.gff3
```
then pulled the sequence of exons using bedtools
```
bedtools getfasta -fi XL9_2.fa -bed XENLA_exon.gff3 -fo exons.fa
```
but this gave me 



