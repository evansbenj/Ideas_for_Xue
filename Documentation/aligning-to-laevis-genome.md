# Genome Alignment

We need to associate the transcripts with the genome at various points.

# Reads direct to the genome

I wanted to explore mapping the trimmed reads directly to the genome, without buidling a transciptome first. The transcriptome build, which will likely have less mapping errors because the transcripts are longer than the reads, will perhaps distinguishing the homeologs better.

In the end, we will likely run differential expression analyses, get values for transcripts, map these back to the genome, and map the transcripts to the genome. For now, I can assess coverage of mapped reads.

*Location:Borealis-Family-Transcriptomes-July2017/Data/Aligned-trimmed-direct-to-genome/*

```bash
# alignment to genome (super scaffold) of each set of reads
../../Scripts/bam-with-readGroups.pl ../Trimmed/ mom # and dad, girl, boy (might split by oviduct and testis)
../../Scripts/bam-with-readGroups.pl ../Trimmed/ boy_testis

# indexing samtools index Mom_sort_$i\.bam Mom_sort_$i\.bai
for i in *dad*bam ; do name=$(grep -o "sort.*\.b" <(echo $i)) ; samtools index $name\am $name\ai ; done
for i in *mom*bam ; do name=$(grep -o "sort.*\.b" <(echo $i)) ; samtools index $name\am $name\ai ; done
for i in *girl_oviduct.bam ; do name=$(grep -o "sort.*\.b" <(echo $i)) ; samtools index $name\am $name\ai ; done
for i in *girl_liver.bam ; do name=$(grep -o "sort.*\.b" <(echo $i)) ; samtools index $name\am $name\ai ; done
for i in *boy_testis.bam ; do name=$(grep -o "sort.*\.b" <(echo $i)) ; samtools index $name\am $name\ai ; done
for i in *boy_liver.bam ; do name=$(grep -o "sort.*\.b" <(echo $i)) ; samtools index $name\am $name\ai ; done


# Genotyping

~/bin/samtools-1.3.1/samtools mpileup -ugf ~/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/Xla.v91.superscaffold.fa -l /home/benf/Borealis_Genome_HiSeqX/Analyses/Mapping_to_laevis/Redo_Sept2016/chr1-2.bed -t DP,AD -s sort*bam | ~/bin/bcftools/bcftools call -V indels --format-fields GQ -m -O z -o mpileup_raw_BorealisTranscriptomes_chr1-2.vcf.gz

~/bin/samtools-1.3.1/samtools mpileup -ugf ~/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/Xla.v91.superscaffold.fa -l /home/benf/Borealis_Genome_HiSeqX/Analyses/Mapping_to_laevis/Redo_Sept2016/chr3-5.bed -t DP,AD -s sort*bam | ~/bin/bcftools/bcftools call -V indels --format-fields GQ -m -O z -o mpileup_raw_BorealisTranscriptomes_chr3-5.vcf.gz

~/bin/samtools-1.3.1/samtools mpileup -ugf ~/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/Xla.v91.superscaffold.fa -l /home/benf/Borealis_Genome_HiSeqX/Analyses/Mapping_to_laevis/Redo_Sept2016/chr6-910.bed -t DP,AD -s sort*bam | ~/bin/bcftools/bcftools call -V indels --format-fields GQ -m -O z -o mpileup_raw_BorealisTranscriptomes_chr6-910.vcf.gz

~/bin/samtools-1.3.1/samtools mpileup -ugf ~/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/Xla.v91.superscaffold.fa -l /home/benf/Borealis_Genome_HiSeqX/Analyses/Mapping_to_laevis/Redo_Sept2016/chr7-8.bed -t DP,AD -s sort*bam | ~/bin/bcftools/bcftools call -V indels --format-fields GQ -m -O z -o mpileup_raw_BorealisTranscriptomes_chr7-8.vcf.gz

```
