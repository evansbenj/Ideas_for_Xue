### Step
1. go through the alignment file and extract only sequence of mapped transcrips
2. redo the quantification with the mapped transcripts
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
kallisto index -i /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/mapped_borealis_transcript_sequence.fa.kallisto_idx /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/mapped_borealis_transcript_sequence.fa

kallisto quant -i /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/mapped_borealis_transcripts/mapped_borealis_transcript_sequence.fa.kallisto_idx  -o female_rep1 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3897_mom_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3897_mom_liver_R2_scythe.fastq.gz)


kallisto quant -i /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/mapped_borealis_transcripts/mapped_borealis_transcript_sequence.fa.kallisto_idx  -o female_rep2 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4009_girl_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4009_girl_liver_R2_scythe.fastq.gz)

kallisto quant -i /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/mapped_borealis_transcripts/mapped_borealis_transcript_sequence.fa.kallisto_idx  -o female_rep3 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4072_girl_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4072_girl_liver_R2_scythe.fastq.gz)

kallisto quant -i /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/mapped_borealis_transcripts/mapped_borealis_transcript_sequence.fa.kallisto_idx  -o female_rep4 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4082_girl_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4082_girl_liver_R2_scythe.fastq.gz)

kallisto quant -i /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/mapped_borealis_transcripts/mapped_borealis_transcript_sequence.fa.kallisto_idx  -o male_rep1 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3896_dad_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3896_dad_liver_R2_scythe.fastq.gz)

kallisto quant -i /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/mapped_borealis_transcripts/mapped_borealis_transcript_sequence.fa.kallisto_idx  -o male_rep2 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3929_boy_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3929_boy_liver_R2_scythe.fastq.gz)

kallisto quant -i /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/mapped_borealis_transcripts/mapped_borealis_transcript_sequence.fa.kallisto_idx  -o male_rep3 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4017_boy_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4017_boy_liver_R2_scythe.fastq.gz)

kallisto quant -i /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/mapped_borealis_transcripts/mapped_borealis_transcript_sequence.fa.kallisto_idx  -o male_rep4 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4039_boy_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4039_boy_liver_R2_scythe.fastq.gz)
```

**collapse the raw count**
```bash

```

**run DE**
```bash

```
