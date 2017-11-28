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
