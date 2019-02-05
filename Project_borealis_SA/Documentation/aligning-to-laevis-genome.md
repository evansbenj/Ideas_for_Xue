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
```

At the end of one of the HiSat runs (dad) gave this message, likely the same for others. Really bad alignment rate. Something to investigate further.

53806825 (100.00%) were paired; of these:
	44448898 (82.61%) aligned concordantly 0 times
	8952875 (16.64%) aligned concordantly exactly 1 time
	405052 (0.75%) aligned concordantly >1 times
	----
	44448898 pairs aligned concordantly 0 times; of these:
		512653 (1.15%) aligned discordantly 1 time
	----
	43936245 pairs aligned 0 times concordantly or discordantly; of these:
		87872490 mates make up the pairs; of these:
			76267187 (86.79%) aligned 0 times
			11047892 (12.57%) aligned exactly 1 time
			557411 (0.63%) aligned >1 times
29.13% overall alignment rate

Checked one of the girls and they also only had a 27% alignment rate. Maybe try tophat?

Can't figure out any obvious answer. One worry was the length of sequence was too short, but that turns out not to be the case as HiSat has no minimum on that. Default is phred33, which is what the reads should be. Maybe I'll try Star for the alignment, but that doesn't answer what it up to produce these poor mapping results. Star wants, but does not require a GTF (like a GFF) but names must match, and my super scaffold reference annotations do not match the names in the GFF (just that one). I could rewrite the GFF to match the super scaffold (anything not in the genome could be dropped from the VCF then just remove the super scaffold...hope it doesn't map...or just go back to the pre-super scaffold reference).

Maybe just try Star to see how effectively it maps the reads. Wish I could figure out HiSat though...

```bash
# Location:/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/
~/bin/STAR-2.5.3a/bin/Linux_x86_64/STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /home/benf/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/STAR-Index/ --genomeFastaFiles /home/benf/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/Xla.v91.superscaffold.fa

# will write into alignment script if works...not liking star though...but see below
~/bin/STAR-2.5.3a/bin/Linux_x86_64/STAR --runThreadN 10 --genomeDir /home/benf/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/STAR-Index/ --readFilesCommand zcat --outFileNamePrefix dad-star --outSAMtype BAM SortedByCoordinate --readFilesIn ../Trimmed/BJE3896_dad_liver_R1_scythe.fastq.gz ../Trimmed/BJE3896_dad_liver_R2_scythe.fastq.gz
```

Had about 54% mapping and 44% not mapping. The not mapping ones were said to be too short, not sure how short that is though. Maybe my trimming was too much? Online suggested not much processing is necessary. So default is 66% matching+penalties of the read length. Maybe I can turn up the number of allowed mismatches for mapping.

Also trying tophat

```bash
bowtie2-build Xla.v91.superscaffold.fa bowtie_index
# should produce a bam file, that I'll need to sort and remove unmapped reads....gave index stupid name I guess, made copy of fa file.
~/bin/tophat-2.1.1.Linux_x86_64/tophat2 -o dad-liver-tophat -p 10   /home/benf/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/bowtie_index ../Trimmed/BJE3896_dad_liver_R1_scythe.fastq.gz ../Trimmed/BJE3896_dad_liver_R2_scythe.fastq.gz
```

Only took a quick look, but apparently only a 2% alignment rate. Going to try and reverse the read order...maybe I am mixing them up.

```bash
~/bin/tophat-2.1.1.Linux_x86_64/tophat2 -o dad-liver-rev-tophat -p 10   /home/benf/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/bowtie_index ../Trimmed/BJE3896_dad_liver_R2_scythe.fastq.gz ../Trimmed/BJE3896_dad_liver_R1_scythe.fastq.gz

~/bin/STAR-2.5.3a/bin/Linux_x86_64/STAR --runThreadN 10 --genomeDir /home/benf/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/STAR-Index/ --readFilesCommand zcat --outFileNamePrefix dad-star-rev --outSAMtype BAM SortedByCoordinate --readFilesIn ../Trimmed/BJE3896_dad_liver_R2_scythe.fastq.gz ../Trimmed/BJE3896_dad_liver_R1_scythe.fastq.gz
```

Reversing the reads had no effect. Star giving the best, try turning down the filters. Will need to check though that it's not substantially worse though.

```bash
# --outFilterMismatchNmax 20 # default is 10
# --outFilterMismatchNoverLmax 0.5 # default is 0.3, ratio mismatches to mapped length
# --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 # length of reads, default is 0.66

~/bin/STAR-2.5.3a/bin/Linux_x86_64/STAR --runThreadN 10 --genomeDir /home/benf/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/STAR-Index/ --readFilesCommand zcat --outFileNamePrefix dad-star-lessStringent --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 20 --outFilterMismatchNoverLmax 0.5 --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 --readFilesIn ../Trimmed/BJE3896_dad_liver_R1_scythe.fastq.gz ../Trimmed/BJE3896_dad_liver_R2_scythe.fastq.gz
```

Increased mapping to 61% and unmapped fell to 34% (too short). Apparently none had too many mismatches.

**Genotyping**
(rerun when read mapping improves.)
```bash
# genotyping...not done yet, has not been necessary...ok running now b/c of sex linked inheritance
~/bin/samtools-1.3.1/samtools mpileup -ugf ~/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/Xla.v91.superscaffold.fa -l /home/benf/Borealis_Genome_HiSeqX/Analyses/Mapping_to_laevis/Redo_Sept2016/chr1-2.bed -t DP,AD -s sort*bam | ~/bin/bcftools/bcftools call -V indels --format-fields GQ -m -O z -o mpileup_raw_BorealisReads_chr1-2.vcf.gz

~/bin/samtools-1.3.1/samtools mpileup -ugf ~/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/Xla.v91.superscaffold.fa -l /home/benf/Borealis_Genome_HiSeqX/Analyses/Mapping_to_laevis/Redo_Sept2016/chr3-5.bed -t DP,AD -s sort*bam | ~/bin/bcftools/bcftools call -V indels --format-fields GQ -m -O z -o mpileup_raw_BorealisReads_chr3-5.vcf.gz

~/bin/samtools-1.3.1/samtools mpileup -ugf ~/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/Xla.v91.superscaffold.fa -l /home/benf/Borealis_Genome_HiSeqX/Analyses/Mapping_to_laevis/Redo_Sept2016/chr6-910.bed -t DP,AD -s sort*bam | ~/bin/bcftools/bcftools call -V indels --format-fields GQ -m -O z -o mpileup_raw_BorealisReads_chr6-910.vcf.gz

~/bin/samtools-1.3.1/samtools mpileup -ugf ~/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/Xla.v91.superscaffold.fa -l /home/benf/Borealis_Genome_HiSeqX/Analyses/Mapping_to_laevis/Redo_Sept2016/chr7-8.bed -t DP,AD -s sort*bam | ~/bin/bcftools/bcftools call -V indels --format-fields GQ -m -O z -o mpileup_raw_BorealisReads_chr7-8.vcf.gz
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

Now call abundances for each. (How is this different than the FPKM and TPKM estimated in step 1? It output files with these numbers.) Note that this is using the merged transcript file and not the gff.

```bash
/home/benf/bin/stringtie-1.3.3b.Linux_x86_64/stringtie -e -G stringtie-merge.out -b mom_liver_ballgown_files -o mom_liver_stringtie_final.gtf ../Aligned-trimmed-direct-to-genome/sort_BJE3897_mom_liver.bam

for i in ../Aligned-trimmed-direct-to-genome/*boy*bam ; do name=$(grep -o "BJE[0-9]*_[a-z]*_[a-z]*" <(echo $i)) ; screen -d -m /home/benf/bin/stringtie-1.3.3b.Linux_x86_64/stringtie -e -G stringtie-merge.out -b $name\_ballgown_files -o $name\_stringtie_final.out $i ; done

for i in ../Aligned-trimmed-direct-to-genome/*girl*bam ; do name=$(grep -o "BJE[0-9]*_[a-z]*_[a-z]*" <(echo $i)) ; screen -d -m /home/benf/bin/stringtie-1.3.3b.Linux_x86_64/stringtie -e -G stringtie-merge.out -b $name\_ballgown_files -o $name\_stringtie_final.out $i ; done

for i in ../Aligned-trimmed-direct-to-genome/*dad*bam ; do name=$(grep -o "BJE[0-9]*_[a-z]*_[a-z]*" <(echo $i)) ; screen -d -m /home/benf/bin/stringtie-1.3.3b.Linux_x86_64/stringtie -e -G stringtie-merge.out -b $name\_ballgown_files -o $name\_stringtie_final.out $i ; done

```

This is all running unbelievably fast, to the point I don't believe it's running right. Or, perhaps it's related to the very low alignment rate seen above. I may have to retry the mapping (this will affect the SNP calling too). I should look at coverage stats, etc., after this (and SNPs) are done. I should also try TopHap (or star) to see if we can get better mapping. May also have to play with HiSat parameters.

Ballgown files have been generated, should check them out, but remember the mapping was terrible.
