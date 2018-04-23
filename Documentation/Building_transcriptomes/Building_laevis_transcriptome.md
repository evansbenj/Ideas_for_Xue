

# Building laevis transcriptome with Trinity (genome-guided method)

### Pre-assembly: mapping reads to laevis genome
Time cost:  92min
```
STAR --runThreadN 20 --genomeDir ~/genome_data/laevis_genome/db_star_laevisGenome_wGTF --readFilesCommand zcat --readFilesIn /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/laevis_liver_R1_scythe.fastq.gz /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/laevis_liver_R2_scythe.fastq.gz --outFileNamePrefix /home/xue/laevis_transcriptome_mar2018/laevis_RNAseq_Star --outSAMtype BAM SortedByCoordinate
```
### Trinity: genome-guided transcriptom assembly
This took about 2.5-3 days. The total number of transcripts is 219207. 
```
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq  --left //home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/laevis_liver_R1_scythe.fastq.gz --right /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/laevis_liver_R2_scythe.fastq.gz --genome_guided_bam /home/xue/laevis_transcriptome_mar2018/laevis_RNAseq_star/laevis_RNAseq_StarAligned.sortedByCoord.out.bam --genome_guided_max_intron 20000 --CPU 25 --full_cleanup --max_memory 200G --output /home/xue/laevis_transcriptome_mar2018/laevis_transcriptome_trinityout; echo "trinity is done at info115 in screen laevis" | mail sarahsongxy@gmail.com
```
### Kallisto
```
/home/xue/software/trinityrnaseq-Trinity-v2.4.0/util/align_and_estimate_abundance.pl --transcripts /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/laevis_genomeguided_transcriptome.fasta --seqType fq --samples_file /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/kallisto_out/laevis_trimmed_samples_files.tsv --est_method kallisto --output_dir /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/kallisto_out --trinity_mode --prep_reference
```
### Covert the quantification into matrix for EdgeR
```
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/util/abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix laevis_gg  --name_sample_by_basedir liver_f1/abundance.tsv liver_f2/abundance.tsv liver_f3/abundance.tsv liver_m1/abundance.tsv liver_m2/abundance.tsv liver_m3/abundance.tsv 
```
### EdgeR: differential analysis 
```
#trinity versin 4.0
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/kallisto_out/laevis_gg.counts.matrix --method edgeR --samples_file /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/kallisto_out/laevis_trimmed_samples_files.tsv --output /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/edgeR_v4_out/

#trinity version 2.2
time /home/xue/software/trinityrnaseq-2.2.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/kallisto_out/laevis_gg.counts.matrix --method edgeR --samples_file /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/kallisto_out/laevis_trimmed_samples_files.tsv --output /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/edgeR_v2_out/

```
### extract the differentially expressed transcripts (~10sec)
Total number of differentially expressed transcripts 
- v4: 
- v2: 453
```
#v4
awk '($4 < -2||$4 >2) && $7<0.05  {print }' /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/edgeR_v4_out/laevis_gg.counts.matrix.female_vs_male.edgeR.DE_results > home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/edgeR_v4_out/laevis_gg_edgeRoutv4_de_filtered.tsv

#v2
awk '($2 < -1||$2 >1) && $5<0.05  {print }' /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/edgeR_v2_out/laevis_gg.counts.matrix.female_vs_male.edgeR.DE_results > /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/edgeR_v2_out/laevis_gg_edgeRoutv2_de_filtered.tsv

```
### extract the sequence of those DE transcripts (~10sec)
```
#v4
perl ~/script/extract_sequence.pl /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/edgeR_v4_out/laevis_gg_edgeRoutv4_de_filtered.tsv /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/laevis_genomeguided_transcriptome.fasta > /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/edgeR_v4_out/laevis_gg_DEtranscript_seq.fasta

#v2
perl ~/script/extract_sequence.pl /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/edgeR_v2_out/laevis_gg_edgeRoutv2_de_filtered.tsv /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/laevis_genomeguided_transcriptome.fasta 1 > /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/edgeR_v2_out/laevis_gg_DEtranscript_seq.fasta
```
### blastn the differential expressed transcripts to the laevis genome

```
#v4
time blastn -task blastn -db ~/genome_data/laevis_genome/db_blastn_laevisGenome/Xl9_2_blastn_db -outfmt 6 -evalue 0.00005 -query /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/edgeR_v4_out/laevis_gg_DEtranscript_seq.fasta -out /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/edgeR_v4_out/DEtranscript_mapping_to_genome_blastout.tsv

#v2
time blastn -task blastn -db ~/genome_data/laevis_genome/db_blastn_laevisGenome/Xl9_2_blastn_db -outfmt 6 -evalue 0.00005 -query /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/edgeR_v2_out/laevis_gg_DEtranscript_seq.fasta -out /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/edgeR_v2_out/DEtranscript_mapping_to_genome_blastout.tsv
```
### filter blastout result and generate a summary
```
#v4
perl ~/script/blastout_filter_summary.pl /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/edgeR_v4_out/DEtranscript_mapping_to_genome_blastout.tsv laevis /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/edgeR_v4_out/tophit_laevis_de_v4

#v2
perl ~/script/blastout_filter_summary.pl /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/edgeR_v2_out/DEtranscript_mapping_to_genome_blastout.tsv laevis /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/edgeR_v2_out/tophit_laevis_de_v2
```
### add the genomic location to the edgeR output of DE transcripts
```
#v4


#v2
perl ~/script/add_location_v3.pl tophit_laevis_de_v2_blastout_tophit.tsv laevis_gg_edgeRoutv2_de_filtered.tsv laevis > laevis_gg_edgeRoutv2_de_filtered_glocation.tsv
```












# Building laevis transcriptome with Trinity (de novo)
### Trinity: *de novo* transcriptome assembly
total number of transcript: 608120
```
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq  --left //home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/laevis_liver_R1_scythe.fastq.gz --right /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/laevis_liver_R2_scythe.fastq.gz --CPU 25 --full_cleanup --max_memory 200G --output /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome_trinityout; echo "trinity is done at info114 in screen laevis" | mail sarahsongxy@gmail.com
```
### Kallisto
```
/home/xue/software/trinityrnaseq-Trinity-v2.4.0/util/align_and_estimate_abundance.pl --transcripts /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/laevis_denovo_transcriptome_trinityout.Trinity.fasta --seqType fq --samples_file /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/laevis_trimmed_samples_files.tsv --est_method kallisto --output_dir /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/kallisto_out --trinity_mode --prep_reference
```
### Covert the quantification into matrix for EdgeR
```
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/util/abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix laevis_denovo  --name_sample_by_basedir liver_f1/abundance.tsv liver_f2/abundance.tsv liver_f3/abundance.tsv liver_m1/abundance.tsv liver_m2/abundance.tsv liver_m3/abundance.tsv 
```
### EdgeR: differential analysis 
```
#trinity versin 4.0
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0//Analysis/DifferentialExpression/run_DE_analysis.pl --matrix /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/kallisto_out/laevis_denovo.counts.matrix --method edgeR --samples_file laevis_trimmed_samples_files.tsv

#trinity version 2.2
time /home/xue/software/trinityrnaseq-2.2.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/kallisto_out/laevis_denovo.counts.matrix --method edgeR --samples_file /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/laevis_trimmed_samples_files.tsv

```
### extract the differentially expressed transcripts (~10sec)
Total number of differentially expressed transcripts 
- v4: 1261
- v2: 1647
```
#v4
awk '($4 < -2||$4 >2) && $7<0.05  {print }' laevis_denovo.counts.matrix.female_vs_male.edgeR.DE_results > laevis_edgeR_result_fdr005.tsv

#v2
awk '($2 < -1||$2 >1) && $5<0.05  {print }' laevis_denovo.counts.matrix.female_vs_male.edgeR.DE_results > laevis_edgeR_result_fdr005.tsv

```
### extract the sequence of those DE transcripts (~10sec)
```
#v4
perl ~/script/extract_sequence.pl /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/edgeR_out/laevis_edgeR_result_fdr005.tsv /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/laevis_denovo_transcriptome_trinityout.Trinity.fasta > laevis_denovo_DEtranscript_seq.fasta

#v2
perl ~/script/extract_sequence.pl /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/edgeR_v2_out/laevis_edgeR_result_fdr005.tsv /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/laevis_denovo_transcriptome_trinityout.Trinity.fasta > /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/edgeR_v2_out/laevis_denovo_DEtranscript_seq.fasta
```

## blastn the differential expressed transcripts to the laevis genome

```
#v4
time blastn -task blastn -db ~/genome_data/laevis_genome/db_blastn_laevisGenome/Xl9_2_blastn_db -outfmt 6 -evalue 0.00005 -query /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/post_DE/laevis_denovo_DEtranscript_seq.fasta -out /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/post_DE/DEtranscript_mapping_to_genome_blastout.tsv

#v2
cd /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/edgeR_v2_out/laevis_denovo_de_blastn/
time blastn -task blastn -db ~/genome_data/laevis_genome/db_blastn_laevisGenome/Xl9_2_blastn_db -outfmt 6 -evalue 0.00005 -query /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/edgeR_v2_out/laevis_denovo_de_blastn/subset_1.fasta -out /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/edgeR_v2_out/laevis_denovo_de_blastn/subset_1_blastout.tsv
```
### filter the blastn output and generate summary
```
```
### combine genomic information and differential expression information
```
cd /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/edgeR_v2_out/laevis_denovo_de_blastn
perl ~/script/add_location_v3.pl laevis_denovo_v2_blastout_tophit.tsv ../laevis_edgeR_result_fdr005.tsv laevis > laevis_edgeR_result_glocation_fromblastn.tsv
```

### mapping those sequence to laevis genome using GMAP (~3min)
```
#v4
gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -A -B 5 -t 8 --cross-species -f samse /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/post_DE/laevis_denovo_DEtranscript_seq.fasta | samtools view -S -b > /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/post_DE/laevis_denovo_DEtranscript_genome_gmap.bam

#v2: 
gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -A -B 5 -t 8 --cross-species -f samse /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/edgeR_v2_out/laevis_denovo_DEtranscript_seq.fasta | samtools view -S -b > /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/edgeR_v2_out/laevis_denovo_DEtranscript_genome_gmap.bam
```
Check on some mapping stat (v2):
- unmapped: 6

### Filter mapping result from GMAP
```
#filter by flag - 0x40 indicated unmapped transcript; and by mapping quality - 0 indicate there is 
samtools view -F 0X4 laevis_denovo_DEtranscript_genome_gmap.bam|awk '$5>0 {print $1,$3,$4}' > laevis_denovo_DEtranscript_mapping_filtered.tsv

#filter
perl filter_gmap_result.pl laevis_denovo_DEtranscript_mapping_filtered.tsv > laevis_denovo_DEtranscript_mapping_filtered2.tsv

#add the genomic location to expression file of DE transcripts
cd /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/edgeR_v2_out
perl add_location_v2.pl laevis_denovo_DEtranscript_mapping_filtered.tsv laevis_edgeR_result_fdr005.tsv > laevis_edgeR_result_fdr005_glocation.tsv
```

### mapping those sequence to borealis DE transcript using Blastn (14sec)
```
#v4
makeblastdb -in /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/edgeR_v2_out/laevis_denovo_DEtranscript_seq.fasta -dbtype nucl -out /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/post_DE/db_laevis_DE_blastn

time blastn -task blastn -db /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/post_DE/db_laevis_DE_blastn -outfmt 6 -evalue 0.00005 -query /home/xue/borealis_DE/liver_mvsf/post_edgeR/liver_trans_fdr005_header.fa -out /home/xue/borealis_DE/liver_mvsf/borealis_laevis_tropicalis_orthologs/borealis_laevis_blastout

#v2
makeblastdb -in /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/edgeR_v2_out/laevis_denovo_DEtranscript_seq.fasta -dbtype nucl -out /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/edgeR_v2_out/db_laevis_DE_blastn

time blastn -task blastn -db /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/edgeR_v2_out/db_laevis_DE_blastn -outfmt 6 -evalue 0.00005 -query /home/xue/borealis_DE/liver_mvsf/post_edgeR/liver_trans_fdr005_header.fa -out /home/xue/borealis_DE/liver_mvsf/borealis_laevis_tropicalis_orthologs/borealis_laevis_v2_blastout.tsv
```

### Star: mapping laevis denovo transcriptome to the laevis genome
```
STARlong --runThreadN 5 --genomeDir /home/xue/genome_data/laevis_genome/db_star_laevisGenome_wGTF --outFileNamePrefix /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/laevis_denovo_TG_mapping_star/laevis_denovo_TG_mapping_star --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 20 --outFilterMismatchNoverLmax 0.5 --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 --outSAMattrRGline ID:Transcriptome SM:allTogether PL:Trinity LB:LB-transcriptome --readFilesIn /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/laevis_denovo_transcriptome_trinityout.Trinity.fasta --seedPerReadNmax 10000

```

