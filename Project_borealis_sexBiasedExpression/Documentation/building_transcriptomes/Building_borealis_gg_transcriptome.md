# Path to files:
```
#path to borealis_genome_v1 from Austin
/home/xue/borealis_transcriptome/borealis_genome_v1/Xbo.v1.fa.gz

#path to borealis RNAseq trimmed data
/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/

#path to borealis genome-guided transcriptome
home/xue/borealis_transcriptome/borealis_gg_transcriptome_may2018/
```

# Trinity: genome-guided transcriptome assembly
```
#Pre-assembly: mapping reads to borealis genome and create a genome_guided_bam file using STAR
#### indexing the borealis genome
STAR --runThreadN 20 --runMode genomeGenerate --genomeDir /home/xue/borealis_transcriptome/borealis_genome_v1/db_borealis_v1_star/ --genomeFastaFiles /home/xue/borealis_transcriptome/borealis_genome_v1/db_borealis_v1_star/Xbo.v1.fa --genomeChrBinNbits 16 
#### mapping borealis RNAseq read to the borealis genome
STAR --runThreadN 20 --runMode alignReads --genomeDir /home/xue/borealis_transcriptome/borealis_genome_v1/db_borealis_v1_star/ --genomeFastaFiles /home/xue/borealis_transcriptome/borealis_genome_v1/Xbo.v1.fa --readFilesCommand zcat --readFilesIn /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3896_dad_liver_R1_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3896_dad_testis_R1_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3897_mom_liver_R1_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3929_boy_liver_R1_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3929_boy_testis_R1_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4009_girl_liver_R1_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4009_girl_oviduct_R1_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4017_boy_liver_R1_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4017_boy_testis_R1_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4039_boy_liver_R1_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4039_boy_testis_R1_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4072_girl_liver_R1_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4072_girl_oviduct_R1_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4082_girl_liver_R1_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4082_girl_oviduct_R1_scythe.fastq.gz /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3896_dad_liver_R2_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3896_dad_testis_R2_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3897_mom_liver_R2_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3929_boy_liver_R2_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE3929_boy_testis_R2_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4009_girl_liver_R2_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4009_girl_oviduct_R2_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4017_boy_liver_R2_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4017_boy_testis_R2_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4039_boy_liver_R2_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4039_boy_testis_R2_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4072_girl_liver_R2_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4072_girl_oviduct_R2_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4082_girl_liver_R2_scythe.fastq.gz,/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4082_girl_oviduct_R2_scythe.fastq.gz --outFileNamePrefix /home/xue/borealis_transcriptome/borealis_genome_v1/borealis_RNAseq_Star --outSAMtype BAM SortedByCoordinate
#### When trying to run the above command, it ran into error message "segmentation fault".  BenF suggest to change --genomeChrBinNbits to 200 and it worked for him before. Star manual suggest suggest --genomeChrBinNbits should be set to  min(18,log2(GenomeLength/NumberOfReferences)). I did the calculation log2(2798753914/23195)~=17 and set --genomeChrBinNbits to 17
STAR --runThreadN 20 --runMode genomeGenerate --genomeDir /home/xue/borealis_transcriptome/borealis_genome_v1/db_borealis_v1_star/ --genomeFastaFiles /home/xue/borealis_transcriptome/borealis_genome_v1/db_borealis_v1_star/Xbo.v1.fa --genomeChrBinNbits 17 
#####I cat all the R1 reads into one file and did the same for R2 reads. I ran Star after the above re-indexing and this one work  
STAR --runThreadN 20 --runMode alignReads --genomeDir /home/xue/borealis_transcriptome/borealis_genome_v1/db_borealis_v1_star/ --readFilesCommand zcat --readFilesIn /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/Borealis_scythe_R1.fastq.gz /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/Borealis_scythe_R2.fastq.gz --outFileNamePrefix /home/xue/borealis_transcriptome/borealis_genome_v1/borealis_RNAseq_Star --outSAMtype BAM SortedByCoordinate
########When I ran the above command, memory max out in the sortedbyCoordinate step (re-run STAR with at least --limitBAMsortRAM 35279761423) so I am going to run it with out sortedbycoordinate and also add a limit on RAM memory usage 
STAR --runThreadN 20 --runMode alignReads --genomeDir /home/xue/borealis_transcriptome/borealis_genome_v1/db_borealis_v1_star/ --readFilesCommand zcat --readFilesIn /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/Borealis_scythe_R1.fastq.gz /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/Borealis_scythe_R2.fastq.gz --outFileNamePrefix /home/xue/borealis_transcriptome/borealis_gg_transcriptome_may2018/borealis_RNAseq_Star --outSAMtype BAM Unsorted

#Trinty assembly (time cost: laevis_ggT_laevis_unigene_blastout_besthit.tsv)
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq  --left /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/Borealis_scythe_R1.fastq.gz --right /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/Borealis_scythe_R2.fastq.gz --genome_guided_bam /home/xue/borealis_transcriptome/borealis_gg_transcriptome_may2018/borealis_RNAseq_StarAligned_coordSorted.bam --genome_guided_max_intron 20000 --CPU 25 --full_cleanup --max_memory 200G --output /home/xue/borealis_transcriptome/borealis_gg_transcriptome_may2018/borealis_gg_transcriptome_trinityout; echo "trinity is done at info114 in screen borealis" | mail sarahsongxy@gmail.com
```
### Kallisto
```bash
/home/xue/software/trinityrnaseq-Trinity-v2.4.0/util/align_and_estimate_abundance.pl --transcripts /home/xue/borealis_transcriptome/borealis_gg_transcriptome_may2018/borealis_gg_transcriptome_trinityout --seqType fq --samples_file /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/kallisto_out/laevis_trimmed_samples_files.tsv --est_method kallisto --output_dir /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/kallisto_out --trinity_mode --prep_reference
```
### Covert the quantification into matrix for EdgeR
```bash
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
- v4: 251
- v2: 453
```
#v4
awk '($4 < -2||$4 >2) && $7<0.05  {print }' /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/edgeR_v4_out/laevis_gg.counts.matrix.female_vs_male.edgeR.DE_results > /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/edgeR_v4_out/laevis_gg_edgeRoutv4_de_filtered.tsv

#v2
awk '($2 < -1||$2 >1) && $5<0.05  {print }' /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/edgeR_v2_out/laevis_gg.counts.matrix.female_vs_male.edgeR.DE_results > /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/edgeR_v2_out/laevis_gg_edgeRoutv2_de_filtered.tsv

```
### extract the sequence of those DE transcripts (~5min)
```
#v4
perl ~/script/extract_sequence.pl /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/edgeR_v4_out/laevis_gg_edgeRoutv4_de_filtered.tsv /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/laevis_genomeguided_transcriptome.fasta 1 > /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/edgeR_v4_out/laevis_gg_DEtranscript_seq.fasta

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
perl ~/script/add_location_v3.pl tophit_laevis_de_v4_blastout_tophit.tsv laevis_gg_edgeRoutv4_de_filtered.tsv laevis > laevis_gg_edgeRoutv4_de_filtered_glocation.tsv

#v2
perl ~/script/add_location_v3.pl tophit_laevis_de_v2_blastout_tophit.tsv laevis_gg_edgeRoutv2_de_filtered.tsv laevis > laevis_gg_edgeRoutv2_de_filtered_glocation.tsv
```
