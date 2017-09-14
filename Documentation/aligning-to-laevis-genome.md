# Genome Alignment

We need to associate the transcripts with the genome at various points.

## Reads direct to the genome

I wanted to explore mapping the trimmed reads directly to the genome, without building a transcriptome first. The transcriptome build, which will likely have less mapping errors because the transcripts are longer than the reads, will perhaps distinguishing the homeologs better.

In the end, we will likely run differential expression analyses based on a *de novo* transcriptome build, get values for transcripts, map these back to the genome, and map the transcripts to the genome. Here though, I will explore read counts on the genome.

###### Need to use splice-aware aligner...but which one?

Splice aware alignments add a specific notation to the cigar for the alignment (N) to indicate that there is a gap in the read alignment. Reading the TOPHAT docs (a popular aligner), they say that TOPHAT has been replaced by HiSAT2, which is supposed to be more accurate and faster. I'll give that a try.

*Location:Borealis-Family-Transcriptomes-July2017/Data/Aligned-trimmed-direct-to-genome/*

```bash
#Location of index: ~/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/HISAT-Index/*
~/bin/hisat2-2.1.0/hisat2-build -p 10 ../Xla.v91.superscaffold.fa hisat-ref

# aligning...apparently should have used the dta option, but only helps with downstream
/home/benf/Borealis-Family-Transcriptomes-July2017/Scripts/align-to-genome-with-readGroups.pl ../Trimmed/ dad --hisat

/home/benf/Borealis-Family-Transcriptomes-July2017/Scripts/align-to-genome-with-readGroups.pl ../Trimmed/ mom --hisat

/home/benf/Borealis-Family-Transcriptomes-July2017/Scripts/align-to-genome-with-readGroups.pl ../Trimmed/ boy_testis --hisat

/home/benf/Borealis-Family-Transcriptomes-July2017/Scripts/align-to-genome-with-readGroups.pl ../Trimmed/ boy_liver --hisat

/home/benf/Borealis-Family-Transcriptomes-July2017/Scripts/align-to-genome-with-readGroups.pl ../Trimmed/ girl_oviduct --hisat

/home/benf/Borealis-Family-Transcriptomes-July2017/Scripts/align-to-genome-with-readGroups.pl ../Trimmed/ girl_liver --hisat

# genotyping...not done yet, has not been necessary
# ~/bin/samtools-1.3.1/samtools mpileup -ugf ~/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/Xla.v91.superscaffold.fa -l /home/benf/Borealis_Genome_HiSeqX/Analyses/Mapping_to_laevis/Redo_Sept2016/chr1-2.bed -t DP,AD -s sort*bam | ~/bin/bcftools/bcftools call -V indels --format-fields GQ -m -O z -o mpileup_raw_BorealisTranscriptomes_chr1-2.vcf.gz
#
# ~/bin/samtools-1.3.1/samtools mpileup -ugf ~/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/Xla.v91.superscaffold.fa -l /home/benf/Borealis_Genome_HiSeqX/Analyses/Mapping_to_laevis/Redo_Sept2016/chr3-5.bed -t DP,AD -s sort*bam | ~/bin/bcftools/bcftools call -V indels --format-fields GQ -m -O z -o mpileup_raw_BorealisTranscriptomes_chr3-5.vcf.gz
#
# ~/bin/samtools-1.3.1/samtools mpileup -ugf ~/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/Xla.v91.superscaffold.fa -l /home/benf/Borealis_Genome_HiSeqX/Analyses/Mapping_to_laevis/Redo_Sept2016/chr6-910.bed -t DP,AD -s sort*bam | ~/bin/bcftools/bcftools call -V indels --format-fields GQ -m -O z -o mpileup_raw_BorealisTranscriptomes_chr6-910.vcf.gz
#
# ~/bin/samtools-1.3.1/samtools mpileup -ugf ~/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/Xla.v91.superscaffold.fa -l /home/benf/Borealis_Genome_HiSeqX/Analyses/Mapping_to_laevis/Redo_Sept2016/chr7-8.bed -t DP,AD -s sort*bam | ~/bin/bcftools/bcftools call -V indels --format-fields GQ -m -O z -o mpileup_raw_BorealisTranscriptomes_chr7-8.vcf.gz
```

## Transcript calling

I am going to follow through with the "new tuxedo" package, which inovlves HiSat (or TopHat), then StringTie (for transcript assembly in the context of a genome), then Ballgown for expression analyses (an R package; or other packages like DESeq).

Need to call transcripts with StringTie.

```bash
/home/benf/bin/stringtie-1.3.3b.Linux_x86_64/stringtie -A mom_liver_gene_abundances.tab -C mom_liver_ref_covered.gtf -G /home/benf/Borealis-Family-Transcriptomes-July2017/Analyses/ZW_Gene_Seqs/XL_9.1_v1.8.3.2.primaryTranscripts.gff3 -o mom_liver_stringtie.out ../Aligned-trimmed-direct-to-genome/sort_BJE3897_mom_liver.bam

for i in ../Aligned-trimmed-direct-to-genome/*boy*bam ; do name=$(grep -o "BJE[0-9]*_[a-z]*_[a-z]*" <(echo $i)) ; screen -d -m  /home/benf/bin/stringtie-1.3.3b.Linux_x86_64/stringtie -A $name\_gene_abundances.tab -C $name\_ref_covered.gtf -G /home/benf/Borealis-Family-Transcriptomes-July2017/Analyses/ZW_Gene_Seqs/XL_9.1_v1.8.3.2.primaryTranscripts.gff3 -o $name\_stringtie.out $i ; done

for i in ../Aligned-trimmed-direct-to-genome/*girl*bam ; do name=$(grep -o "BJE[0-9]*_[a-z]*_[a-z]*" <(echo $i)) ; screen -d -m  /home/benf/bin/stringtie-1.3.3b.Linux_x86_64/stringtie -A $name\_gene_abundances.tab -C $name\_ref_covered.gtf -G /home/benf/Borealis-Family-Transcriptomes-July2017/Analyses/ZW_Gene_Seqs/XL_9.1_v1.8.3.2.primaryTranscripts.gff3 -o $name\_stringtie.out $i ; done

for i in ../Aligned-trimmed-direct-to-genome/*dad*bam ; do name=$(grep -o "BJE[0-9]*_[a-z]*_[a-z]*" <(echo $i)) ; screen -d -m  /home/benf/bin/stringtie-1.3.3b.Linux_x86_64/stringtie -A $name\_gene_abundances.tab -C $name\_ref_covered.gtf -G /home/benf/Borealis-Family-Transcriptomes-July2017/Analyses/ZW_Gene_Seqs/XL_9.1_v1.8.3.2.primaryTranscripts.gff3 -o $name\_stringtie.out $i ; done
```
That step doesn't seem right, as it took minutes to run! But, sequences from all chrs and scaffolds seem to be there.

Merge transcripts
```bash
~/bin/stringtie-1.3.3b.Linux_x86_64/stringtie --merge -G /home/benf/Borealis-Family-Transcriptomes-July2017/Analyses/ZW_Gene_Seqs/XL_9.1_v1.8.3.2.primaryTranscripts.gff3 -o stringtie-merge.out *out
```
