# Evolution of expression directionality
### Purpose
The question we can to answer:
- does DE transcript/genes in the sex-link region 


# Finding Orthologs
### general method 
1. bidirectional best hit/reciprocal best hit
  - blast the list of borealis DE (ex, transcript B) into laevis transcriptome => ex, hit transcript L
  - blast transcript L into borealis transcriptome
  - if L's best hit is B => ortholog; if not => not ortholog
2. identify orthology by genomic region
  - GMAP/STAR borealis DE into laevis genome -> done; path ''
  - GMAP/STAR laevis transcriptome into laevis genome
  - orthology = borealis transcript that is in the same genomic region of laevis transcript

# Reciprocal best hit
### borealis-laevis orthologs - all DE transcript
finding ortholog through reciprocal best hit
```
#blast borealis DE to laevis transcriptome (time cost: 6min)
time blastn -task blastn -db /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/db_laevis_denovo_transcriptome -outfmt 6 -evalue 0.00005 -query /home/xue/borealis_DE/liver_mvsf/post_edgeR/borealis_liver_de_transcriptSeq.fa -out /home/xue/borealis_DE/liver_mvsf/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs/borealisDE_laevisTranscriptome_blastout.tsv

#identify the tophit of each transcript; tophit (lowest evalue, highest bitscore)
perl ~/script/blastout_filter_summary.pl borealisDE_laevisTranscriptome_blastout.tsv laevis tophit_borealisDE_laevisT
rm -f *problematic*
rm -f *summary*

# extract seq of the tophit laevis transcripts 
perl ~/script/extract_sequence.pl /home/xue/borealis_DE/liver_mvsf/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs/borealisDE_laevisT_blastout_tophit.tsv /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/laevis_denovo_transcriptome_trinityout.Trinity.fasta 2 > /home/xue/borealis_DE/liver_mvsf/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs/bl_tophit_laevis_transcriptSeq.fa

#building blastn db index for borealis transcriptome
makeblastdb -in /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta -dbtype nucl -out /home/xue/borealis_DE/db_borealist_transcriptome_blastn/db_borealist_transcriptome_blastn

#blast the tophit laevis seq to borealis transcriptome
time blastn -task blastn -db /home/xue/borealis_DE/db_borealis_transcriptome_blastn/db_borealis_transcriptome_blastn -outfmt 6 -evalue 0.00005 -query /home/xue/borealis_DE/liver_mvsf/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs/bl_tophit_laevis_transcriptSeq.fa -out /home/xue/borealis_DE/liver_mvsf/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs/laevisTophit_borealiseT_blastout.tsv

#identify the tophit of each transcript; tophit (lowest evalue, highest bitscore)
cd /home/xue/borealis_DE/liver_mvsf/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs/
perl ~/script/blastout_filter_summary.pl laevisTophit_borealiseT_blastout.tsv laevis tophit_laevisTophit_borealisT
rm -f *summary*


#with a perl script: check if transcripts are reciprocal best hit for each other -> need a perl script to do this
perl ~/script/ortholog_check_rbh.pl borealisDE_laevisT_blastout_tophit.tsv laevisTophit_borealisT_blastout_tophit.tsv > ortholog_borealis_de_laevis_denovoT_list.tsv
```
Extract expression level in TPM for each transcript
```
#extract 
```



Just out of curiousity: if gmap and blastn will get the same result. Theoritically they should be very similar
GMAP: mapping borealis DE to laevis_denovo_transcriptome; 
```
#indexing the laevis_denovo_transcriptome
gmap_build -D /home/xue/laevis_transcriptome_mar2018/db_laevis_denovo_transcriptome_gmap -s none -d laevis_denovo_transcriptome_gmap   /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/laevis_denovo_transcriptome_trinityout.Trinity.fasta

#indexing the laevis_gg_transcriptome
gmap_build -D /home/xue/laevis_transcriptome_mar2018/db_laevis_gg_transcriptome_gmap -s none -d laevis_gg_transcriptome_gmap   /home/xue/laevis_transcriptome_mar2018/laevis_gg_transcriptome/laevis_genomeguided_transcriptome.fasta

#mapping borealis de transcripts to laevis_denovo_transcriptome: 96 unmapped (unmapped defined by 0x4 flag in field #2); 
gmap -D /home/xue/laevis_transcriptome_mar2018/db_laevis_denovo_transcriptome_gmap -d laevis_denovo_transcriptome_gmap  -A -B 5 -t 8 -f samse --cross-species /home/xue/borealis_DE/liver_mvsf/post_edgeR/borealis_liver_de_transcriptSeq.fa | samtools view -S -b > /home/xue/borealis_DE/liver_mvsf/borealis_laevis_tropicalis_orthologs/borealis_de_laevis_transcriptome_gmap/borealis_de_laevis_denovoT_gmap.bam


#mapping borealis de transcripts to laevis_gg_transcriptome: 278 mapped (unmapped defined by 0x4 flag in field #2); 
gmap -D /home/xue/laevis_transcriptome_mar2018/db_laevis_gg_transcriptome_gmap -d laevis_gg_transcriptome_gmap -A -B 5 -t 8 -f samse --cross-species /home/xue/borealis_DE/liver_mvsf/post_edgeR/borealis_liver_de_transcriptSeq.fa | samtools view -S -b > /home/xue/borealis_DE/liver_mvsf/borealis_laevis_tropicalis_orthologs/borealis_de_laevis_transcriptome_gmap/borealis_de_laevis_ggT_gmap.bam
```


# Identify orthologs by genomic location
### borealis-laevis orthologs
blastn: mapping borealis DE and laevis_gg_transcriptom to laevis genome 
 - blast is good for looking for highly similar sequence but might map the most conserve region of the transcript to the most converse region of the gene and might not identify the right gene
```
cd /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/blastn_transcriptome_genome/subset_seq
perl ~/script/split_fasta.pl ../../laevis_genomeguided_transcriptome.fasta 9000

time blastn -task blastn -db ~/genome_data/laevis_genome/db_blastn_laevisGenome/Xl9_2_blastn_db -outfmt 6 -evalue 0.00005 -query /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/blastn_transcriptome_genome/subset_seq/subset_25.fasta -out /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/blastn_transcriptome_genome/subset_blastout/subset_25_blastout.tsv
```
STARlong: mapping borealis DE and laevis_gg_transcriptom to laevis genome
```
#mapping borealis_de_transcripts to laevis genome: only 29 transcripts were mapped; it is weird, what if I try without all the filtering parameter 
STARlong --runThreadN 5 --genomeDir /home/xue/genome_data/laevis_genome/db_star_laevisGenome_wGTF --outFileNamePrefix /home/xue/borealis_DE/liver_mvsf/mapping_star/borealis_de_laevis_genome_mapping_star --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 20 --outFilterMismatchNoverLmax 0.5 --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 --outSAMattrRGline ID:Transcriptome SM:allTogether PL:Trinity LB:LB-transcriptome --readFilesIn /home/xue/borealis_DE/liver_mvsf/mapping_star/liver_DEseq_IDonly.fasta --seedPerReadNmax 10000

#mapping borealis_de_transcripts to laevis genome: without all the filtering parameter, nothing mapped; re-run the above
STARlong --runThreadN 5 --genomeDir /home/xue/genome_data/laevis_genome/db_star_laevisGenome_wGTF --outFileNamePrefix /home/xue/borealis_DE/liver_mvsf/mapping_star/borealis_de_laevis_genome_mapping_star --outSAMtype BAM SortedByCoordinate  --readFilesIn /home/xue/borealis_DE/liver_mvsf/mapping_star/liver_DEseq_IDonly.fasta --seedPerReadNmax 10000

#mapping laevis_gg_transcriptom to laevis genome
STARlong --runThreadN 5 --genomeDir /home/xue/genome_data/laevis_genome/db_star_laevisGenome_wGTF --outFileNamePrefix /home/xue/borealis_DE/liver_mvsf/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_byGenomicLocation/laevis_gg_TG_mapping_star --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 20 --outFilterMismatchNoverLmax 0.5 --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 --outSAMattrRGline ID:Transcriptome SM:allTogether PL:Trinity LB:LB-transcriptome --readFilesIn /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/laevis_genomeguided_transcriptome.fasta --seedPerReadNmax 10000
```
GMAP: mapping borealis DE and laevis_denovo_transcriptom to laevis genome
```
#mapping borealis_de_transcripts to laevis genome:already done;
gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -A -Z -f samse /home/xue/borealis_DE/liver_mvsf/filtered_edgeRout/liver_trans_fdr005_header.fa > /home/xue/borealis_DE/liver_mvsf/mapping_GMAP/liver_DE_gmap_out.sam
#the sam file were converted to bam file and store at:
/home/xue/borealis_DE/liver_mvsf/mapping_GMAP/liver_DE_gmap_out.bam
cd /home/xue/borealis_DE/liver_mvsf/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_byGenomicLocation
mv /home/xue/borealis_DE/liver_mvsf/mapping_GMAP/liver_DE_gmap_out.bam . 


#mapping laevis_denovo_transcriptom to laevis genome (time cost: ~15hrs)
gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -A -B 5 -t 8 -f samse /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/laevis_denovo_transcriptome_trinityout.Trinity.fasta | samtools view -S -b > /home/xue/borealis_DE/liver_mvsf/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_byGenomicLocation/mapped_gmap/laevis_denovo_transcriptome_genome_gmap.bam
```
Bedtool: extract alignment coordinates for transcripts based on their CIGAR strings; 
```
bedtools bamtobed -i liver_DE_gmap_out.bam > liver_DE_gmap_out_bedfile.tsv
bedtools bamtobed -i laevis_denovo_transcriptome_genome_gmap.bam > laevis_denovo_transcriptome_genome_gmap_bedfile.tsv
```
Identify orthology by running perl script orthologs_identification_byGenomicLocation.pl
```
perl ~/script/orthologs_identification_byGenomicLocation.pl liver_DE_gmap_out_bedfile.tsv laevis_denovo_transcriptome_genome_gmap_bedfile.tsv > list_orthologs_borealis_de_laevis_transcriptome.tsv
```

GMAP: mapping borealis DE and laevis_denovo_transcriptome to laevis v91 genome; just because the laevis-tropicalis alignment files contains alignment to laevis v91 genome
```
#the folder for the output are in 
/home/xue/borealis_DE/liver_mvsf/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_byGenomicLocation/mapped_gmap_xla_v91

#indexing the laevis_v91 genome
gmap_build -D /home/xue/genome_data/laevis_genome/db_gmap_xl91 -s none -g -d laevis91_gmap ../Xla.v91.fa.gz

#mapping borealis_de_transcripts to laevis_v91 genome
gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl91 -d laevis91_gmap -A -B 5 -t 8 -f samse --cross-species /home/xue/borealis_DE/liver_mvsf/post_edgeR/borealis_liver_de_transcriptSeq.fa | samtools view -S -b > /home/xue/borealis_DE/liver_mvsf/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_byGenomicLocation/mapped_gmap_xla_v91/borealis_de_laevis91_genome_gmap.bam

#mapping laevis_denovo_transcriptome to laevis_v91 genome
gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl91 -d laevis91_gmap -A -B 5 -t 8 -f samse /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/laevis_denovo_transcriptome_trinityout.Trinity.fasta | samtools view -S -b > /home/xue/borealis_DE/liver_mvsf/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_byGenomicLocation/mapped_gmap_xla_v91/laevis_denovo_transcriptome_genome_gmap.bam
```


### laevis-tropicalis orthologs
We got the alignment files from Austin.
The files are stored @
 `/home/benf/Borealis-Family-Transcriptomes-July2017/Data/genomic-alignment-tropicalis-laevisLS`
 and a back-up copy @
 `/home/xue/borealis_DE/genomic-alignment-tropicalis-laevisLS`

The formatting of the file is:
- TotalT.tab is a table file where the first column is the alignment "name" (starts at aln_0 and counts up), second is the "class" where XLS is a 3-way alignment between trop, laevis L, and laevis S, XL is just trop and L, XS trop and S, LS is only L and S. The rest of the columns are the scaffold, start, stop strand of the sequences based on the class (meaning class LS lists the laevis L location first, and has no information for tropicalis)

GMAP: mapping tropicalis transcriptome to tropicalis v71 genome; just because the laevis-tropicalis alignment files contains alignment of tropicalis v71 genome and laevis_v91 genome
```
#indexing the tropicalis_v71_genome
gmap_build -D /home/xue/genome_data/tropicalis_genome/db_gmap_tropicalisv71  -s none -g -d tropicalis71_gmap /home/xue/genome_data/tropicalis_genome/tropicalis_v71_main_genome/Xenopus_tropicalis_main_genome_scaffolds.fasta.gz

#mapping tropicali_denovo_transcriptome to tropicalis_v71_genome
gmap -D /home/xue/genome_data/tropicalis_genome/db_gmap_tropicalisv71 -d tropicalis71_gmap -A -B 5 -t 8 -f samse /home/xue/tropicalis_transcriptome/tropicalis_denovo_transcriptome/tropicalis_trinityout.Trinity.fasta | samtools view -S -b > /home/xue/tropicalis_transcriptome/tropicalis_denovo_transcriptome/tropicalis_denovoT_genomev71_mapping_gmap/tropicalis_denovo_transcriptome_genomeV71_gmap.bam


#mapping tropicali_gg_transcriptome to tropicalis_v71_genome
gmap -D /home/xue/genome_data/tropicalis_genome/db_gmap_tropicalisv71 -d tropicalis71_gmap -A -B 5 -t 8 -f samse /home/xue/tropicalis_transcriptome/tropicalis_gg_transcriptome/tropicalis_gg_transcriptome.fasta | samtools view -S -b > /home/xue/tropicalis_transcriptome/tropicalis_gg_transcriptome/tropicalis_ggT_genomev71_mapping_gmap/tropicalis_gg_transcriptome_genomeV71_gmap.bam

```

### borealis-laevis-tropicalis orthologs

# Identifying orthologs: combination method
The method will combine partial output from bbh and genomic-based method. 
Basically, I will align borealis de transcripts to laevis denovo transcriptome. For the best hits, look at their genomic location. The genomic location of borealis de transcripts is obtained by mapping de transcripts to laevis genome v92. The genomic locaiton of laevis transcripts is obtained by mapping laevis transcripts to laevis genome v92. Most of those information are available from previous work.

#### path to files
borealis de transcripts to laevis denovo transcriptome: 
  - blastn
  `/home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_bbh/borealisDE_laevisTranscriptome_blastout.tsv`
  - gmap
  `/home/xue/borealis_DE/liver_mvsf/borealis_laevis_tropicalis_orthologs/borealis_de_laevis_transcriptome_gmap/borealis_de_laevis_denovoT_gmap.bam`

borealis de transcript to laevis genome v92 
```
/home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_byGenomicLocation/mapped_gmap/liver_DE_gmap_out_bedfile.tsv
```

laevis denovo transcriptome to laevis genome v92
  ```
  /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_byGenomicLocation/mapped_gmap/laevis_denovo_transcriptome_genome_gmap_bedfile.tsv
  ```

#### perl script
```
#working directory
cd /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_byGenomicLocation/mapped_gmap

#select blast best hit for each DE transcripts
perl ~/script/blastout_besthit.pl /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_bbh/borealisDE_laevisTranscriptome_blastout.tsv > borealis_de_laevis_denovoT_blastout_besthit.tsv 

#
perl ~/script/orthologs_TT_TG_output_analysis.pl borealis_de_laevis_denovoT_blastout_besthit.tsv ../laevis_denovo_transcriptome_genome_gmap_bedfile.tsv ../liver_DE_gmap_out_bedfile.tsv|less

```


