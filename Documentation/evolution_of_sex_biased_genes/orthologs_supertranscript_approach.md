# Orthologs identification

directory:`/home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach`

## Check supertranscript
split and map then check if the parts mapped to same chromosome
```
cd /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach
perl ~/script/split_seq.pl SuperDuper.fasta > SuperDuper_splitted.fasta 

```
### mapping to borealis genome
mapping: gmap
```
gmap -D /home/xue/borealis_transcriptome/borealis_genome_v1/db_borealis_genome_v1_gmap -d borealis_genome_v1 -A -B 5 -t 8 -f samse /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/SuperDuper_splitted.fasta | samtools view -S -b > /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/SuperDuper_splitted_gmap.bam
```

mapping: star
```
STARlong --runThreadN 5 --genomeDir /home/xue/borealis_transcriptome/borealis_genome_v1/db_borealis_v1_star --outFileNamePrefix /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/SuperDuper_splitted_star --outSAMtype BAM Unsorted  --readFilesIn /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/SuperDuper_splitted.fasta --seedPerReadNmax 10000
```
Checking the mapping result and count the number of problematic transcripts, which are transcripts that have parts mapping to different chromosomes.
```
cd /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/mapped_to_borealis_genome/

perl ~/script/check_supertranscript.pl SuperDuper_splitted_gmap.bam > SuperDuper_splitted_gmap_summary_noLS.tsv
vi ~/script/check_supertranscript.pl #command in the chop to chop off the 'L' or 'S'
perl ~/script/check_supertranscript.pl SuperDuper_splitted_gmap.bam > SuperDuper_splitted_gmap_summary_LS.tsv
```
### mapping to laevis unigene
mapping: gmap
```
# making gmap database for laevis unigene
gmap_build -d db_gmap_laevisUnigene/ -D /home/xue/genome_data/laevis_genome/ -g ../Xl.seq.uniq.gz 

gmap -d db_gmap_laevisUnigene/ -D /home/xue/genome_data/laevis_genome/ -A -B 5 -t 8 -f samse --cross-species /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/SuperDuper_splitted.fasta | samtools view -S -b > /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/SuperDuper_splitted_gmap_laevisUnigene.bam
```

mapping: star
```
#indexing laevis unigene for star mapping
STAR --runThreadN 20 --runMode genomeGenerate --genomeDir /home/xue/genome_data/laevis_genome/db_star_laevisUnigene --genomeFastaFiles /home/xue/genome_data/laevis_genome/Xl_seq_uniq.fa --genomeChrBinNbits 16

#mapping splitted supertanscript to laevis unigene
STARlong --runThreadN 5 --genomeDir /home/xue/genome_data/laevis_genome/db_star_laevisUnigene --outFileNamePrefix /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/mapped_to_laevis_unigene/SuperDuper_splitted_star_laevisUnigene --outSAMtype BAM Unsorted --readFilesIn /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/SuperDuper_splitted.fasta --seedPerReadNmax 10000
```
Checking the mapping result and count the number of problematic transcripts, which are transcripts that have parts mapping to different chromosomes.
```
cd /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/mapped_to_laevis_unigene

vi ~/script/check_supertranscript.pl #command out the if check for 'hr'
perl ~/script/check_supertranscript.pl SuperDuper_splitted_gmap_laevisUnigene.bam > SuperDuper_splitted_gmap_laevisUnigene_summary_LS.tsv
perl ~/script/check_supertranscript.pl SuperDuper_splitted_star_laevisUnigeneAligned.out.bam > SuperDuper_splitted_star_laevisUnigene_summary_LS.tsv
vi ~/script/check_supertranscript.pl #command in the chop to chop off the 'L' or 'S'
perl ~/script/check_supertranscript.pl SuperDuper_splitted_gmap_laevisUnigene.bam > SuperDuper_splitted_gmap_laevisUnigene_summary_noLS.tsv
perl ~/script/check_supertranscript.pl SuperDuper_splitted_star_laevisUnigeneAligned.out.bam > SuperDuper_splitted_star_laevisUnigene_summary_noLS.tsv
```
### mapping to laevis genome
mapping: gmap
```
gmap -d laevis92_gmap/ -D /home/xue/genome_data/laevis_genome/db_gmap_xl92 -A -B 5 -t 8 -f samse --cross-species /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/SuperDuper_splitted.fasta | samtools view -S -b > /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/mapped_to_laevis_genome/SuperDuper_splitted_gmap_laevisGenome.bam
```

mapping: star
```
#mapping splitted supertanscript to laevis genome
STARlong --runThreadN 5 --genomeDir /home/xue/genome_data/laevis_genome/db_star_laevisGenome_wGTF --outFileNamePrefix /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/mapped_to_laevis_genome/SuperDuper_splitted_star_laevisGenome --outSAMtype BAM Unsorted --readFilesIn /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/SuperDuper_splitted.fasta --seedPerReadNmax 10000
```
Checking the mapping result and count the number of problematic transcripts, which are transcripts that have parts mapping to different chromosomes.
```
cd /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/mapped_to_laevis_genome/

perl ~/script/check_supertranscript.pl SuperDuper_splitted_gmap_laevisGenome.bam > SuperDuper_splitted_gmap_laevisGenome_summary_LS.tsv
perl ~/script/check_supertranscript.pl SuperDuper_splitted_star_laevisGenomeAligned.out.bam > SuperDuper_splitted_star_laevisGenome_summary_LS.tsv
vi ~/script/check_supertranscript.pl #command in the chop to chop off the 'L' or 'S'
perl ~/script/check_supertranscript.pl SuperDuper_splitted_gmap_laevisGenome.bam > SuperDuper_splitted_gmap_laevisGenome_summary_noLS.tsv
perl ~/script/check_supertranscript.pl SuperDuper_splitted_star_laevisGenomeAligned.out.bam > SuperDuper_splitted_star_laevisGenome_summary_noLS.tsv
```

## borealis_laevis orthologs
Use borealis supertranscript as ref genome. Map everything to it
  - indexing the supertranscripts
  ```
  #indexing for blastn
  makeblastdb -in /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/SuperDuper.fasta  -dbtype nucl -out /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/db_blastn_borealisSuper/blastn_borealisSuper

  #indexing for gmap
  gmap_build -d db_gmap_borealisSuper/ -D /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/SuperDuper.fasta 

  #indexing for star
  STAR --runThreadN 20 --runMode genomeGenerate --genomeDir /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/db_star_borealisSuper --genomeFastaFiles /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/SuperDuper.fasta --genomeChrBinNbits 16 --limitGenomeGenerateRAM=64660815232
  ```
  - mapping: borealis de transcript to borealis super-transcripts 
  ```
  #blastn borealis_De to borealis_supertranscript
  time blastn -task blastn -db /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/db_blastn_borealisSuper/blastn_borealisSuper -outfmt 6 -evalue 0.00005 -query /home/xue/borealis_DE/de_sex_liver/post_edgeR/borealis_liver_de_transcriptSeq.fa -out /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/supertranscript_as_refGenome/borealis_de_borealis_superT_blastout.tsv

  #gmap borealis_De to borealis_supertranscript
  gmap -d db_gmap_borealisSuper/ -D /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach -A -Z -f samse /home/xue/borealis_DE/de_sex_liver/post_edgeR/borealis_liver_de_transcriptSeq.fa | samtools view -S -b >  /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/supertranscript_as_refGenome/borealis_de_borealis_superT_gmap.bam

  ```
  - mapping: laevis transcripts to borealis super-transcripts 
  ```
  #blastn laevis ggT to borealis_supertranscript
  time blastn -task blastn -db /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/db_blastn_borealisSuper/blastn_borealisSuper -outfmt 6 -evalue 0.00005 -query /home/xue/laevis_transcriptome_mar2018/laevis_gg_transcriptome/laevis_genomeguided_transcriptome.fasta -out /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/supertranscript_as_refGenome/laevis_ggT_borealis_superT_blastout.tsv

  #blastn laevis ggT to borealis_supertranscript
  time blastn -task blastn -db /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/db_blastn_borealisSuper/blastn_borealisSuper -outfmt 6 -evalue 0.00005 -query /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/laevis_denovo_transcriptome_trinityout.Trinity.fasta -out /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/supertranscript_as_refGenome/laevis_denovoT_borealis_superT_blastout.tsv

  #gmap laevis ggT to borealis_supertranscript
  gmap -d db_gmap_borealisSuper/ -D /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach -A -Z -f samse --cross-species /home/xue/laevis_transcriptome_mar2018/laevis_gg_transcriptome/laevis_genomeguided_transcriptome.fasta | samtools view -S -b >  /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/supertranscript_as_refGenome/laevis_ggT_borealis_superT_gmap.bam


  #gmap laevis denovo transcriptome to borealis_supertranscript
  gmap -d db_gmap_borealisSuper/ -D /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach -A -Z -f samse --cross-species /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/laevis_denovo_transcriptome_trinityout.Trinity.fasta | samtools view -S -b >  /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_supertranscriptApproach/supertranscript_as_refGenome/laevis_denovoT_borealis_superT_gmap.bam
  ```
   - mapping: tropicalis transcripts to borealis super-transcripts 
  ```
 #blastn tropicalis ggT to borealis_supertranscript - screen gg
  time blastn -task blastn -db /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_supertranscriptApproach/db_blastn_borealisSuper/blastn_borealisSuper -outfmt 6 -evalue 0.00005 -query /home/xue/tropicalis_transcriptome/tropicalis_gg_transcriptome/tropicalis_gg_transcriptome.fasta -out /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_supertranscriptApproach/tropicalis_ggT_borealis_superT_blastout.tsv

  #blastn tropicalis denovoT to borealis_supertranscript - screen denovo
  time blastn -task blastn -db /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_supertranscriptApproach/db_blastn_borealisSuper/blastn_borealisSuper -outfmt 6 -evalue 0.00005 -query /home/xue/tropicalis_transcriptome/tropicalis_denovo_transcriptome/tropicalis_denovo_transcriptome_trinityout.fasta -out /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_supertranscriptApproach/supertranscript_as_refGenome/tropicalis_denovoT_borealis_superT_blastout.tsv

  #gmap tropicalis ggT to borealis_supertranscript - screen tropicalis
  gmap -d db_gmap_borealisSuper/ -D /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_supertranscriptApproach -A -Z -f samse --cross-species /home/xue/tropicalis_transcriptome/tropicalis_gg_transcriptome/tropicalis_gg_transcriptome.fasta | samtools view -S -b >  /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_supertranscriptApproach/supertranscript_as_refGenome/tropicalis_ggT_borealis_superT_gmap.bam


  #gmap tropicalis denovo transcriptome to borealis_supertranscript - screen laevis
  gmap -d db_gmap_borealisSuper/ -D /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_supertranscriptApproach -A -Z -f samse --cross-species /home/xue/tropicalis_transcriptome/tropicalis_denovo_transcriptome/tropicalis_denovo_transcriptome_trinityout.fasta | samtools view -S -b > /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_supertranscriptApproach/supertranscript_as_refGenome/tropicalis_denovoT_borealis_superT_gmap.bam
  ```

Filter the mapping output
  -  filter: keep top hit only (lowest evalue, highest bitscore)
  ```
  #borealis de
  #perl ~/script/blastout_tophit.pl borealis_de_borealis_superT_blastout.tsv > borealis_de_borealis_superT_blastout_tophit.tsv
  perl ~/script/blastout_tophit.pl borealis_de_borealis_superT_blastout.tsv borealis_de_borealis_superT
  #laevis ggT
  #perl ~/script/blastout_tophit.pl laevis_ggT_borealis_superT_blastout.tsv > laevis_ggT_borealis_superT_blastout_tophit.tsv
  perl ~/script/blastout_tophit.pl laevis_ggT_borealis_superT_blastout.tsv laevis_ggT_borealis_superT
  ```
  - Combine the two tophit files
  ```
  #add the "xb_" and "xl_" in front of the trinity transcript id of borealis and laevis respectively. 
  sed 's/^/xb_/' borealis_de_borealis_superT_blastout_tophit.tsv > borealis_de_borealis_superT_blastout_tophit_flagged.tsv;
  sed 's/^/xl_/' laevis_ggT_borealis_superT_blastout_tophit.tsv > laevis_ggT_borealis_superT_blastout_tophit_flagged.tsv

  #combine them together
  cat borealis_de_borealis_superT_blastout_tophit_flagged.tsv laevis_ggT_borealis_superT_blastout_tophit_flagged.tsv > borealis_laevis_transcriptome_borealis_superT_tophit_flagged.tsv
  ```
  - sort the file and group the transcript IDs by unigene IDs
  ```
  sort -k 2,2 borealis_laevis_transcriptome_borealis_superT_tophit_flagged.tsv > borealis_laevis_transcriptome_borealis_superT_tophit_flagged_sorted.tsv
  ```
  - keep only transcript that mapped to supertranscripts that match a borealis de transcript
  ``` 
   perl ~/script/xl_xb_to_unigene_blastout_filter.pl borealis_laevis_transcriptome_borealis_superT_tophit_flagged_sorted.tsv > borealis_laevis_transcriptome_borealis_superT_tophit_flagged_sorted_filtered.tsv
  
  ```
