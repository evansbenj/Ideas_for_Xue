# Identification of orthologs using laevis and tropicalis unigene
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
blastn: blast borealis DE and laevis transcriptome into laevis unigenes
```
#blast borealis DE into laevis unigenes
time blastn -task blastn -db /home/xue/genome_data/laevis_genome/db_blastn_laevisUnigene/db_blastn_laevisUnigene -outfmt 6 -evalue 0.00005 -query /home/xue/borealis_DE/de_sex_liver/post_edgeR/borealis_liver_de_transcriptSeq.fa -out /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/borealis_laevis_orthologs_unigeneApproach/borealis_de_laevis_unigene_blastout.tsv 

#blast laevis transcriptome into laevis unigenes
cd /home/xue/laevis_transcriptome_mar2018/laevis_gg_transcriptome

time blastn -task blastn -db /home/xue/genome_data/laevis_genome/db_blastn_laevisUnigene/db_blastn_laevisUnigene -outfmt 6 -evalue 0.00005 -query /home/xue/laevis_transcriptome_mar2018/laevis_gg_transcriptome/subset/subset_8.fasta -out /home/xue/laevis_transcriptome_mar2018/laevis_gg_transcriptome/subset/blastout_subset_8.fasta
```
