### Steps
1. go through the alignment file and extract only sequence of mapped transcrips *done
2. redo the quantification with the mapped transcripts *running
3. collapse the raw count
4. run DE

**working directory**
```bash
/home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/mapped_borealis_transcripts
```

**extract sequence of the mapped transcripts**
```bash
perl ~/script/extract_sequence.pl borealis_denovoT_laevisV92_genome_gmap.bed trinity_out_dir.Trinity.fasta 4 >mapped_borealis_transcript_sequence.fa
```

**quantification**
```bash
#Indexing the transcriptomes
kallisto index -i /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/mapped_borealis_transcript_sequence.fa.kallisto_idx /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/mapped_borealis_transcript_sequence.fa

#Transcript abundance quantification
kallisto quant -i /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/mapped_borealis_transcripts/mapped_borealis_transcript_sequence.fa.kallisto_idx  -o female_rep1 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3897_mom_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3897_mom_liver_R2_scythe.fastq.gz)


kallisto quant -i /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/mapped_borealis_transcripts/mapped_borealis_transcript_sequence.fa.kallisto_idx  -o female_rep2 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4009_girl_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4009_girl_liver_R2_scythe.fastq.gz)

kallisto quant -i /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/mapped_borealis_transcripts/mapped_borealis_transcript_sequence.fa.kallisto_idx  -o female_rep3 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4072_girl_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4072_girl_liver_R2_scythe.fastq.gz)

kallisto quant -i /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/mapped_borealis_transcripts/mapped_borealis_transcript_sequence.fa.kallisto_idx  -o female_rep4 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4082_girl_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4082_girl_liver_R2_scythe.fastq.gz)

kallisto quant -i /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/mapped_borealis_transcripts/mapped_borealis_transcript_sequence.fa.kallisto_idx  -o male_rep1 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3896_dad_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3896_dad_liver_R2_scythe.fastq.gz)

kallisto quant -i /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/mapped_borealis_transcripts/mapped_borealis_transcript_sequence.fa.kallisto_idx  -o male_rep2 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3929_boy_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3929_boy_liver_R2_scythe.fastq.gz)

kallisto quant -i /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/mapped_borealis_transcripts/mapped_borealis_transcript_sequence.fa.kallisto_idx  -o male_rep3 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4017_boy_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4017_boy_liver_R2_scythe.fastq.gz)

kallisto quant -i /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/mapped_borealis_transcripts/mapped_borealis_transcript_sequence.fa.kallisto_idx  -o male_rep4 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4039_boy_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4039_boy_liver_R2_scythe.fastq.gz)

#To compute the matrix (time cost: ~10min):
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/util/abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix borealis_liver  --name_sample_by_basedir female_rep1/abundance.tsv female_rep2/abundance.tsv female_rep3/abundance.tsv female_rep4/abundance.tsv male_rep1/abundance.tsv male_rep2/abundance.tsv male_rep3/abundance.tsv male_rep4/abundance.tsv

```

**collapse the raw count**
```bash
time perl ~/script/identify_genomic_location_for_transcripts.pl --matrix /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/mapped_borealis_transcripts/borealis_liver.counts.matrix > borealis_expression_matrix_per_gene.tsv
```

**run DE**
```bash

```


## Assess the quality of the reference transcriptome
```
#build a mapping index
bowtie2-build ../mapped_borealis_transcript_sequence.fa ../mapped_borealis_transcript_sequence.fa 

#
bowtie2 -p 10 -q --no-unal -k 20 -x Trinity.fasta -1 /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/Borealis_scythe_R1.fastq.gz -2 /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/Borealis_scythe_R2.fastq.gz align_stats.txt| samtools view -@10 -Sb -o reads_alignment_to_mapped_borealis_transcript_bowtie2.bam 
     
cat 2>&1 align_stats.txt

```

# laevis
**working directory**

```bash
cd /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/collapsed_laevis_transriptome_DE/countQuantification_mappedTransOnly
```

**extract sequence of the mapped transcripts**
```bash
perl ~/script/extract_sequence.pl /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/laevis_denovo_TG_mapping_gmap/laevis_denovoT_laevisV92_genome_gmap_bedfile_filtered.bed /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/laevis_denovo_transcriptome_trinityout.Trinity.fasta 4 > laevis_denovo_transcriptome_mapped_sequence.fa
```

**quantification**
```bash
#Indexing the transcriptomes
kallisto index -i /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/collapse_laevis_transcriptome_withMappedTransOnly/laevis_denovo_transcriptome_mapped_sequence.fa.kallisto_idx /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/collapse_laevis_transcriptome_withMappedTransOnly/laevis_denovo_transcriptome_mapped_sequence.fa

#Transcript abundance quantification
kallisto quant -i /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/collapse_laevis_transcriptome_withMappedTransOnly/laevis_denovo_transcriptome_mapped_sequence.fa.kallisto_idx  -o female_rep1 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/BJE4524_girl_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/BJE4524_girl_liver_R2_scythe.fastq.gz)


kallisto quant -i /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/collapse_laevis_transcriptome_withMappedTransOnly/laevis_denovo_transcriptome_mapped_sequence.fa.kallisto_idx  -o female_rep2 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/BJE4525_girl_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/BJE4525_girl_liver_R2_scythe.fastq.gz)

kallisto quant -i /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/collapse_laevis_transcriptome_withMappedTransOnly/laevis_denovo_transcriptome_mapped_sequence.fa.kallisto_idx  -o female_rep3 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/BJE4526_girl_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/BJE4526_girl_liver_R2_scythe.fastq.gz)

kallisto quant -i /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/collapse_laevis_transcriptome_withMappedTransOnly/laevis_denovo_transcriptome_mapped_sequence.fa.kallisto_idx  -o male_rep1 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/BJE4527_boy_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/BJE4527_boy_liver_R2_scythe.fastq.gz)

kallisto quant -i /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/collapse_laevis_transcriptome_withMappedTransOnly/laevis_denovo_transcriptome_mapped_sequence.fa.kallisto_idx  -o male_rep2 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/BJE4528_boy_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/BJE4528_boy_liver_R2_scythe.fastq.gz)

kallisto quant -i /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/collapse_laevis_transcriptome_withMappedTransOnly/laevis_denovo_transcriptome_mapped_sequence.fa.kallisto_idx  -o male_rep3 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/BJE4529_boy_liver_R1_scythe.fastq.gz ) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/BJE4529_boy_liver_R2_scythe.fastq.gz )


#To compute the matrix (time cost: ~10min):
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/util/abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix laevis_liver  --name_sample_by_basedir female_rep1/abundance.tsv female_rep2/abundance.tsv female_rep3/abundance.tsv  male_rep1/abundance.tsv male_rep2/abundance.tsv male_rep3/abundance.tsv 

```

**collapse the raw count**
```bash
time perl ~/script/identify_genomic_location_for_transcripts.pl --matrix /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/mapped_borealis_transcripts/borealis_liver.counts.matrix > borealis_expression_matrix_per_gene.tsv
```



