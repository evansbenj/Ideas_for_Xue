

# Building laevis transcriptome with Trinity (genome-guided method)

### Pre-assembly: mapping reads to laevis genome
Time cost:  92min
```
STAR --runThreadN 20 --genomeDir ~/genome_data/laevis_genome/db_star_laevisGenome_wGTF --readFilesCommand zcat --readFilesIn /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/laevis_liver_R1_scythe.fastq.gz /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/laevis_liver_R2_scythe.fastq.gz --outFileNamePrefix /home/xue/laevis_transcriptome_mar2018/laevis_RNAseq_Star --outSAMtype BAM SortedByCoordinate
```
### Trinity: genome-guided transcriptom assembly

```
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq  --left //home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/laevis_liver_R1_scythe.fastq.gz --right /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Trimmed/laevis_liver_R2_scythe.fastq.gz --genome_guided_bam /home/xue/laevis_transcriptome_mar2018/laevis_RNAseq_star/laevis_RNAseq_StarAligned.sortedByCoord.out.bam --genome_guided_max_intron 20000 --CPU 25 --full_cleanup --max_memory 200G --output /home/xue/laevis_transcriptome_mar2018/laevis_transcriptome_trinityout; echo "trinity is done at info115 in screen laevis" | mail sarahsongxy@gmail.com
```
# Building laevis transcriptome with Trinity (de novo)
### Trinity: *de novo* transcriptome assembly
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
### differential analysis with EdgeR
```
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0//Analysis/DifferentialExpression/run_DE_analysis.pl --matrix /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/kallisto_out/laevis_denovo.counts.matrix --method edgeR --samples_file laevis_trimmed_samples_files.tsv
```
### extract the differential expressed transcripts
Total number: 11806
```
awk '($2 < -1||$2 >1) && $5<0.05  {print }' laevis_denovo.counts.matrix.female_vs_male.edgeR.DE_results > laevis_fdr005.tsv

```
### extract the sequence of those DE transcripts
```
perl ~/script/extract_sequence.pl /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/edgeR_out/laevis_fdr005.tsv /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/laevis_denovo_transcriptome_trinityout.Trinity.fasta > laevis_denovo_DEtranscript_seq.fasta
```
### mapping those sequence to laevis genome using GMAP
```
gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -A -Z -f samse /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/post_DE/laevis_denovo_DEtranscript_seq.fasta > /home/xue/borealis_DE/liver_mvsf/mapping_GMAP/liver_DE_gmap_out.sam

gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -A -S -B 5 -t 8 --cross-species -S /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/post_DE/laevis_denovo_DEtranscript_seq.fasta > alignment_only

gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -A -B 5 -t 8 --cross-species -f samse /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/post_DE/laevis_denovo_DEtranscript_seq.fasta | samtools view -S -b > /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/post_DE/laevis_denovo_DEtranscript_genome_gmap.bam
```

### mapping those sequence to borealis DE transcript using Blastn


