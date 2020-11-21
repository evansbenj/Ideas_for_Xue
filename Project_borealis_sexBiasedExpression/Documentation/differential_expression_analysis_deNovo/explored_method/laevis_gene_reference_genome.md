I would like to build a laevis reference genome, which each gene contained sequence of its exons only. If there are two exons that are overlapping (splice variants), the longer exons will be kept. 
```bash
#extract exon coordinate information from gff file
awk '$3 == "exon"{print }' XENLA_9.2_Xenbase.gff3 > XENLA_exon.gff3

#extract the sequence of exons
 bedtools getfasta -fi ../XL9_2.fa -bed XENLA_exon.gff3 -fo XENLA_exon.fa
 
 #use gffread to extract the sequence of gene
gffread/gffread/gffread XENLA_9.2_Xenbase.gff3 -g XL9_2.fa  -ZM -w xl_mRNA.fasta 
 
```
