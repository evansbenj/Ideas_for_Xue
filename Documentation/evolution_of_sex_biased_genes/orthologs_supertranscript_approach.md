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


