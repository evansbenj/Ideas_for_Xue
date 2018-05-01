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
GMAP: mapping borealis DE and laevis_gg_transcriptom to laevis genome
```
#mapping borealis_de_transcripts to laevis genome:already done; the path is /home/xue/borealis_DE/liver_mvsf/mapping_GMAP/liver_DE_gmap_out.bam
cd /home/xue/borealis_DE/liver_mvsf/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_byGenomicLocation
mv /home/xue/borealis_DE/liver_mvsf/mapping_GMAP/liver_DE_gmap_out.bam . 


#mapping laevis_gg_transcriptom to laevis genome
gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -A -Z -f samse /home/xue/borealis_DE/liver_mvsf/filtered_edgeRout/liver_trans_fdr005_header.fa > /home/xue/borealis_DE/liver_mvsf/mapping_GMAP/liver_DE_gmap_out.sam
```

### laevis-tropicalis orthologs
We got the alignment files from Austin.
The files are stored in two different account:
 `/home/benf/Borealis-Family-Transcriptomes-July2017/Data/genomic-alignment-tropicalis-laevisLS`
 and 
 `/home/xue/borealis_DE/genomic-alignment-tropicalis-laevisLS`

The formatting of the file is:
    - TotalT.tab is a table file where the first column is the alignment "name" (starts at aln_0 and counts up), second is the "class" where XLS is a 3-way alignment between trop, laevis L, and laevis S, XL is just trop and L, XS trop and S, LS is only L and S. The rest of the columns are the scaffold, start, stop strand of the sequences based on the class (meaning class LS lists the laevis L location first, and has no information for tropicalis)
```
```


### borealis-laevis-tropicalis orthologs


