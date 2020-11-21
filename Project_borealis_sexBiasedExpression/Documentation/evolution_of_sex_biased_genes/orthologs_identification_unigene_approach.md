# Identification of orthologs using laevis and tropicalis unigene
working directory:
```
/home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_unigeneApproach
```
getting unigene 
```
#getting tropicalis unigenes
cd /home/xue/genome_data/tropicalis_genome
wget ftp://ftp.ncbi.nih.gov/repository/UniGene/Xenopus_tropicalis/Str.seq.uniq.gz

#getting laevis unigenes
cd /home/xue/genome_data/laevis_genome
wget ftp://ftp.ncbi.nih.gov/repository/UniGene/Xenopus_laevis/Xl.seq.uniq.gz
```
Blastn indexing of unigenes
```
#indexing tropicalis unigenes
makeblastdb -in /home/xue/genome_data/tropicalis_genome/str_seq_uniq.fa -dbtype nucl -out /home/xue/genome_data/tropicalis_genome/db_blastn_tropicalisUnigene/db_blastn_tropicalisUnigene

#indexing laevis unigenes
makeblastdb -in /home/xue/genome_data/laevis_genome/Xl_seq_uniq.fa -dbtype nucl -out /home/xue/genome_data/laevis_genome/db_blastn_laevisUnigene/db_blastn_laevisUnigene
```
## borealis_laevis orthologs
blastn: blast borealis DE and laevis transcriptome into laevis unigenes
```
#blast borealis DE into laevis unigenes
time blastn -task blastn -db /home/xue/genome_data/laevis_genome/db_blastn_laevisUnigene/db_blastn_laevisUnigene -outfmt 6 -evalue 0.00005 -query /home/xue/borealis_DE/de_sex_liver/post_edgeR/borealis_liver_de_transcriptSeq.fa -out /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_unigeneApproach/borealis_de_laevis_unigene_blastout.tsv 

#blast laevis transcriptome into laevis unigenes
cd /home/xue/laevis_transcriptome_mar2018/laevis_gg_transcriptome

time blastn -task blastn -db /home/xue/genome_data/laevis_genome/db_blastn_laevisUnigene/db_blastn_laevisUnigene -outfmt 6 -evalue 0.00005 -query /home/xue/laevis_transcriptome_mar2018/laevis_gg_transcriptome/subset/subset_8.fasta -out /home/xue/laevis_transcriptome_mar2018/laevis_gg_transcriptome/laevis_ggT_laevis_unigene_blastout_besthit.tsv
```
filter: keep best hit only (lowest evalue)
```
#borealis de
perl ~/script/blastout_besthit.pl borealis_de_laevis_unigene_blastout.tsv > borealis_de_laevis_unigene_blastout_besthit.tsv

#laevis transcriptome
scp /home/xue/laevis_transcriptome_mar2018/laevis_gg_transcriptome/laevis_ggT_laevis_unigene_blastn/laevis_ggT_laevis_unigene_blastout_besthit.tsv .

```
Combine the two besthit files
```
#add the "xb_" and "xl_" in front of the trinity transcript id of borealis and laevis respectively. 
sed 's/^/xb_/' borealis_de_laevis_unigene_blastout_besthit.tsv > borealis_de_laevis_unigene_blastout_besthit_flagged.tsv
sed 's/^/xl_/' laevis_ggT_laevis_unigene_blastout_besthit.tsv > laevis_ggT_laevis_unigene_blastout_besthit_flagged.tsv

#combine them together
cat borealis_de_laevis_unigene_blastout_besthit.tsv laevis_ggT_laevis_unigene_blastout_besthit.tsv > borealis_laevis_transcriptome_laevis_unigene_besthit.tsv

cat borealis_de_laevis_unigene_blastout_besthit_flagged.tsv laevis_ggT_laevis_unigene_blastout_besthit_flagged.tsv > borealis_laevis_transcriptome_laevis_unigene_besthit_flagged.tsv
```
sort the file and group the transcript IDs by unigene IDs
```
sort -k 2,2 borealis_laevis_transcriptome_laevis_unigene_besthit.tsv > borealis_laevis_transcriptome_laevis_unigene_besthit_sorted_filtered.tsv
sort -k 2,2 borealis_laevis_transcriptome_laevis_unigene_besthit_flagged.tsv > borealis_laevis_transcriptome_laevis_unigene_besthit_flagged_sorted_filtered.tsv
```
extracting the sequence of unique laevis unigenes
```
sort -u -k2,2 borealis_laevis_transcriptome_laevis_unigene_besthit_sorted_filtered.tsv |awk '$2 {print}' > borealis_laevis_transcriptome_laevis_unigene_besthit_unigeneID.tsv

perl ~/script/extract_sequence.pl borealis_laevis_transcriptome_laevis_unigene_besthit_unigeneID.tsv ~/genome_data/laevis_genome/Xl_seq_uniq.fa 1 > borealis_laevis_transcriptome_laevis_unigene_besthit_unigeneSeq.fa
```
gmap the laevis unigenes to the laevis genome to find out their genomic location
```
gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -A -B 5 -t 8 -f samse /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_unigeneApproach/borealis_laevis_transcriptome_laevis_unigene_besthit_unigeneSeq.fa | samtools view -S -b > /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_unigeneApproach/borealis_laevis_transcriptome_laevis_unigene_besthit_unigeneSeq_gmap.bam

```
After sorting the catanated files, one borealis transcript and multiple laevis transcripts might map to one laevis unigene transcript. So a thing to try is to map borealis_DE transcript to laevis_genomeguided_transcriptome and then pick only the top hit pair. 
```
#indexing laevis_ggT for blastn mapping;
makeblastdb -in /home/xue/laevis_transcriptome_mar2018/laevis_gg_transcriptome/laevis_genomeguided_transcriptome.fasta -dbtype nucl -out /home/xue/laevis_transcriptome_mar2018/db_laevis_gg_transcriptome_blastn/db_blastn_laevis_ggT

time blastn -task blastn -db /home/xue/laevis_transcriptome_mar2018/db_laevis_gg_transcriptome_blastn/db_blastn_laevis_ggT -outfmt 6 -evalue 0.00005 -query /home/xue/borealis_DE/de_sex_liver/post_edgeR/borealis_liver_de_transcriptSeq.fa -out /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_unigeneApproach/borealis_de_laevis_ggT_blastout.tsv 


```
Another thing to try is to pick the laevis transcript that mapped to laevis unigene with the lowest e-value and highest bitscore
```
#current label position: laevis_ggT_ID  laevis_unigene_ID 

#pick the top hit for each unigene ID

```
out of curiousity, does two xb_de transcript mapped to the same unigene? -> no 
```
sort -u -k2,2 borealis_de_laevis_unigene_blastout_besthit_flagged.tsv > borealis_de_laevis_unigene_blastout_besthit_flagged_sorted.tsv

awk 'x[$2]++ == 1 { print $2 " is duplicated"}' borealis_de_laevis_unigene_blastout_besthit_flagged_sorted.tsv
```

## laevis_tropicalis orthologs
blast laevis unigene and tropicalis unigenes to each other
```
#blast laevis to tropicalis
blastn -task megablast -db /home/xue/genome_data/tropicalis_genome/db_blastn_tropicalisUnigene/db_blastn_tropicalisUnigene -query /home/xue/genome_data/laevis_genome/Xl_seq_uniq.fa -outfmt 6 -evalue 1e-20 -max_target_seqs 2 -out /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_unigeneApproach/ST_to_XL_unigene_blastnout.tsv 

#blast tropicalis to laevis
blastn -task megablast -db /home/xue/genome_data/laevis_genome/db_blastn_laevisUnigene/db_blastn_laevisUnigene -query /home/xue/genome_data/tropicalis_genome/str_seq_uniq.fa -outfmt 6 -evalue 1e-20 -max_target_seqs 2 -out /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_unigeneApproach/ST_to_XL_unigene_megablastout.tsv

#blast tropicalis transcriptome to tropicalis unigene
blastn -task megablast -db /home/xue/genome_data/tropicalis_genome/db_blastn_tropicalisUnigene/db_blastn_tropicalisUnigene -query /home/xue/tropicalis_transcriptome/tropicalis_gg_transcriptome/tropicalis_gg_transcriptome.fasta -outfmt 6 -evalue 1e-20 -max_target_seqs 2 -out /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_unigeneApproach/st_transcriptome_to_st_unigene_blastnout.tsv

blastn -task megablast -db /home/xue/genome_data/tropicalis_genome/db_blastn_tropicalisUnigene/db_blastn_tropicalisUnigene -query /home/xue/tropicalis_transcriptome/tropicalis_denovo_transcriptome/tropicalis_trinityout.Trinity.fasta -outfmt 6 -evalue 1e-20 -max_target_seqs 2 -out /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_unigeneApproach/st_denovoT_to_st_unigene_blastnout.tsv
```
run Makes_diads_and_triads_XLunigene_STensembl.pl from BenE. This generate output ST_unigene_XL_unigene_dyads_and_triads_unigene_minlen_200 and trop_unigene_IDs_with_no_reciprocal_best_blast_hit_unigene
```
cd /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_unigeneApproach
perl ~/script/Makes_diads_and_triads_XLunigene_STensembl.pl XL_to_ST_unigene_megablastout.tsv ST_to_XL_unigene_megablastout.tsv
```
identify tropicalis unigene orthologs of laevis unigenes that is orthology to borealis DE
```
cd /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_unigeneApproach
perl ~/script/orthologs_st_xl_trias_dias_filter.pl borealis_laevis_transcriptome_laevis_unigene_besthit_unigeneID.tsv ST_unigene_XL_unigene_dyads_and_triads_unigene_minlen_200 >  ST_unigene_XL_unigene_dyads_and_triads_borealisDE_ortho.tsv
```
filter tropicalis blastn out besthit 
```
#only select transcripts that mapped to tropicalis unigene orthologs
perl ~/script/orthologs_st_xl_trias_dias_filter.pl ST_unigene_XL_unigene_dyads_and_triads_borealisDE_ortho.tsv st_ggT_to_st_unigene_blastnout_besthit.tsv > st_ggT_to_st_unigene_blastnout_besthit_borealisDE_ortho.tsv
#only select keep the tophit
perl ~/script/blastout_tophit.pl st_ggT_to_st_unigene_blastnout_besthit_borealisDE_ortho.tsv  st_ggT_to_st_unigene_blastnout_besthit_borealisDE_ortho
```
filter tropicalis blastn out besthit 
```
#only select transcripts that mapped to laevis unigene orthologs
perl ~/script/orthologs_st_xl_trias_dias_filter.pl borealis_laevis_transcriptome_laevis_unigene_besthit_unigeneID.tsv laevis_ggT_laevis_unigene_blastout_besthit.tsv > laevis_ggT_laevis_unigene_blastout_besthit_borealis_ortho.tsv
#only select keep the tophit
perl ~/script/blastout_tophit.pl laevis_ggT_laevis_unigene_blastout_besthit_borealis_ortho.tsv laevis_ggT_laevis_unigene_blastout_besthit_borealis_ortho

```
## XB_XL_ST orthologs
Use a perl script to compile the ortholog information together and extract the ID of transcripts. Logic:
```
#idenntifying orthology - unigene approach
#combining all the mapping output and extract the ortholog information

#read file borealis_laevis_transcriptome_laevis_unigene_besthit.tsv
#store the information in a hash 
#$xb_xl_ortho{XB_DE_ID}{XL_unigene}{laevis_t_ID/st_unigene_ID/ST_T_ID}[] = ID_name
#
```
