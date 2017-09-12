# Genome Alignment

We need to associate the transcripts with the genome at various points.

# Reads direct to the genome

I wanted to explore mapping the trimmed reads directly to the genome, without building a transcriptome first. The transcriptome build, which will likely have less mapping errors because the transcripts are longer than the reads, will perhaps distinguishing the homeologs better.

In the end, we will likely run differential expression analyses based on a *de novo* transcriptome build, get values for transcripts, map these back to the genome, and map the transcripts to the genome. Here though, I will explore read counts on the genome.

###### Need to use splice-aware aligner...but which one?

Splice aware alignments add a specific notation to the cigar for the alignment (N) to indicate that there is a gap in the read alignment. Reading the TOPHAT docs (a popular aligner), they say that TOPHAT has been replaced by HiSAT2, which is supposed to be more accurate and faster. I'll give that a try.

*Location:Borealis-Family-Transcriptomes-July2017/Data/Aligned-trimmed-direct-to-genome/*

```bash
#Location of index: ~/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/HISAT-Index/*
~/bin/hisat2-2.1.0/hisat2-build -p 10 ../Xla.v91.superscaffold.fa hisat-ref

# aligning
/home/benf/Borealis-Family-Transcriptomes-July2017/Scripts/align-to-genome-with-readGroups.pl ../Trimmed/ dad --hisat

/home/benf/Borealis-Family-Transcriptomes-July2017/Scripts/align-to-genome-with-readGroups.pl ../Trimmed/ mom --hisat

/home/benf/Borealis-Family-Transcriptomes-July2017/Scripts/align-to-genome-with-readGroups.pl ../Trimmed/ boy_testis --hisat

/home/benf/Borealis-Family-Transcriptomes-July2017/Scripts/align-to-genome-with-readGroups.pl ../Trimmed/ boy_liver --hisat

/home/benf/Borealis-Family-Transcriptomes-July2017/Scripts/align-to-genome-with-readGroups.pl ../Trimmed/ girl_oviduct --hisat

/home/benf/Borealis-Family-Transcriptomes-July2017/Scripts/align-to-genome-with-readGroups.pl ../Trimmed/ girl_liver --hisat

# genotyping
~/bin/samtools-1.3.1/samtools mpileup -ugf ~/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/Xla.v91.superscaffold.fa -l /home/benf/Borealis_Genome_HiSeqX/Analyses/Mapping_to_laevis/Redo_Sept2016/chr1-2.bed -t DP,AD -s sort*bam | ~/bin/bcftools/bcftools call -V indels --format-fields GQ -m -O z -o mpileup_raw_BorealisTranscriptomes_chr1-2.vcf.gz

~/bin/samtools-1.3.1/samtools mpileup -ugf ~/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/Xla.v91.superscaffold.fa -l /home/benf/Borealis_Genome_HiSeqX/Analyses/Mapping_to_laevis/Redo_Sept2016/chr3-5.bed -t DP,AD -s sort*bam | ~/bin/bcftools/bcftools call -V indels --format-fields GQ -m -O z -o mpileup_raw_BorealisTranscriptomes_chr3-5.vcf.gz

~/bin/samtools-1.3.1/samtools mpileup -ugf ~/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/Xla.v91.superscaffold.fa -l /home/benf/Borealis_Genome_HiSeqX/Analyses/Mapping_to_laevis/Redo_Sept2016/chr6-910.bed -t DP,AD -s sort*bam | ~/bin/bcftools/bcftools call -V indels --format-fields GQ -m -O z -o mpileup_raw_BorealisTranscriptomes_chr6-910.vcf.gz

~/bin/samtools-1.3.1/samtools mpileup -ugf ~/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/Xla.v91.superscaffold.fa -l /home/benf/Borealis_Genome_HiSeqX/Analyses/Mapping_to_laevis/Redo_Sept2016/chr7-8.bed -t DP,AD -s sort*bam | ~/bin/bcftools/bcftools call -V indels --format-fields GQ -m -O z -o mpileup_raw_BorealisTranscriptomes_chr7-8.vcf.gz
```
