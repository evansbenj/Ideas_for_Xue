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

### extract exons using bedtool
Ideally, I would like to pull transcripts and each of them contains exons of the gene. It seems like bedtools doesn't have the option to group the exons by gene and then pull the gene sequence. So i did some initial filtering and create a gff file that contains only the exons.
```
grep -w "exon" XENLA_9.2_Xenbase.gff3 > XENLA_exon.gff3
```
then pulled the sequence of exons using bedtools
```
bedtools getfasta -fi XL9_2.fa -bed XENLA_exon.gff3 -fo exons.fa
```
but this gave me 

