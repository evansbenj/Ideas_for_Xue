# DE analysis between sex using all tissue samples
### Kallisto
```
time perl /home/xue/software/trinityrnaseq-Trinity-v2.4.0/util/align_and_estimate_abundance.pl --transcripts /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta --seqType fa --samples_file borealis_de_all_samplefile.tsv --est_method kallisto --output_dir /home/xue/borealis_DE/all_de_sex/borealis_all_kallisto_denovo --trinity_mode --prep_reference
```
### Covert the quantification into matrix for EdgeR (~10min)
```
cd /home/xue/borealis_DE/all_de_sex/borealis_all_kallisto_denovo
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/util/abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix /home/xue/borealis_DE/all_de_sex/borealis_all_kallisto_denovo/borealis_all  --name_sample_by_basedir female_rep1/abundance.tsv female_rep2/abundance.tsv female_rep3/abundance.tsv female_rep4/abundance.tsv female_rep5/abundance.tsv female_rep6/abundance.tsv female_rep7/abundance.tsv male_rep1/abundance.tsv male_rep2/abundance.tsv male_rep3/abundance.tsv male_rep4/abundance.tsv male_rep5/abundance.tsv male_rep6/abundance.tsv male_rep7/abundance.tsv
```
### EdgeR: differential analysis 
```
#trinity versin 4.0
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix /home/xue/borealis_DE/all_de_sex/borealis_all_kallisto_denovo/borealis_all.counts.matrix --method edgeR --samples_file /home/xue/borealis_DE/all_de_sex/borealis_all_kallisto_denovo/borealis_de_all_samplefile.tsv --output /home/xue/borealis_DE/all_de_sex/borealis_all_edgeR_v4_denovo/

#trinity version 2.2 (~13min)
time /home/xue/software/trinityrnaseq-2.2.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix /home/xue/borealis_DE/all_de_sex/borealis_all_kallisto_denovo/borealis_all.counts.matrix --method edgeR --samples_file /home/xue/borealis_DE/all_de_sex/borealis_all_kallisto_denovo/borealis_de_all_samplefile.tsv --output /home/xue/borealis_DE/all_de_sex/borealis_all_edgeR_v2_denovo/

```
### extract the differentially expressed transcripts (~10sec)
Total number of differentially expressed transcripts 
- v4: 
- v2: 961
```
#v4
awk '($4 < -2||$4 >2) && $7<0.05  {print }' /home/xue/borealis_DE/all_de_sex/borealis_all_edgeR_v4_denovo/borealis_all.counts.matrix.female_vs_male.edgeR.DE_results > /home/xue/borealis_DE/all_de_sex/borealis_all_edgeR_v4_denovo/borealis_all_edgeRoutv4_de_filtered.tsv

#v2 (~3sec)
awk '($2 < -1||$2 >1) && $5<0.05  {print }' /home/xue/borealis_DE/all_de_sex/borealis_all_edgeR_v2_denovo/borealis_all.counts.matrix.female_vs_male.edgeR.DE_results > /home/xue/borealis_DE/all_de_sex/borealis_all_edgeR_v2_denovo/borealis_all_edgeRoutv2_de_filtered.tsv

```
### extract the sequence of those DE transcripts 
```
#v4
perl ~/script/extract_sequence.pl /home/xue/borealis_DE/all_de_sex/borealis_all_edgeR_v4_denovo/borealis_all_edgeRoutv4_de_filtered.tsv /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta 1 > /home/xue/borealis_DE/all_de_sex/borealis_all_edgeR_v4_denovo/borealis_all_DEtranscript_seq.fasta

#v2 (~1.15min)
perl ~/script/extract_sequence.pl /home/xue/borealis_DE/all_de_sex/borealis_all_edgeR_v2_denovo/borealis_all_edgeRoutv2_de_filtered.tsv /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta 1 > /home/xue/borealis_DE/all_de_sex/borealis_all_edgeR_v2_denovo/borealis_all_DEtranscript_seq.fasta
```
### blastn the differential expressed transcripts to the laevis genome 

```
#v4
time blastn -task blastn -db ~/genome_data/laevis_genome/db_blastn_laevisGenome/Xl9_2_blastn_db -outfmt 6 -evalue 0.00005 -query /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/edgeR_v4_out/laevis_gg_DEtranscript_seq.fasta -out /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/edgeR_v4_out/DEtranscript_mapping_to_genome_blastout.tsv

#v2 (~2hrs)
time blastn -task blastn -db ~/genome_data/laevis_genome/db_blastn_laevisGenome/Xl9_2_blastn_db -outfmt 6 -evalue 0.00005 -query /home/xue/borealis_DE/all_de_sex/borealis_all_edgeR_v2_denovo/borealis_all_DEtranscript_seq.fasta -out /home/xue/borealis_DE/all_de_sex/borealis_all_edgeR_v2_denovo/DEtranscript_mapping_to_genome_blastout.tsv
```
### filter blastout result and generate a summary
```
#v4
perl ~/script/blastout_filter_summary.pl /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/edgeR_v4_out/DEtranscript_mapping_to_genome_blastout.tsv laevis /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/edgeR_v4_out/tophit_laevis_de_v4

#v2
perl ~/script/blastout_filter_summary.pl /home/xue/borealis_DE/all_de_sex/borealis_all_edgeR_v2_denovo/DEtranscript_mapping_to_genome_blastout.tsv borealis /home/xue/borealis_DE/all_de_sex/borealis_all_edgeR_v2_denovo/borealis_all_de_v2
```
### add the genomic location to the edgeR output of DE transcripts
```
#v4
perl ~/script/add_location_v3.pl tophit_laevis_de_v4_blastout_tophit.tsv laevis_gg_edgeRoutv4_de_filtered.tsv laevis > laevis_gg_edgeRoutv4_de_filtered_glocation.tsv

#v2
perl ~/script/add_location_v3.pl tophit_laevis_de_v2_blastout_tophit.tsv laevis_gg_edgeRoutv2_de_filtered.tsv laevis > laevis_gg_edgeRoutv2_de_filtered_glocation.tsv
```
