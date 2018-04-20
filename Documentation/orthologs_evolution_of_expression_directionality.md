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

```
Extract expression level in TPM for each transcript
```
#extract 
```




  
# Identify orthologs by genomic location
### borealis-laevis orthologs
blastn: mapping laevis transcriptome to laevis genome; 
 - blast is good for looking for highly similar sequence but might map the most conserve region of the transcript to the most converse region of the gene and might not identify the right gene
```
cd /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/blastn_transcriptome_genome/subset_seq
perl ~/script/split_fasta.pl ../../laevis_genomeguided_transcriptome.fasta 9000

time blastn -task blastn -db ~/genome_data/laevis_genome/db_blastn_laevisGenome/Xl9_2_blastn_db -outfmt 6 -evalue 0.00005 -query /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/blastn_transcriptome_genome/subset_seq/subset_25.fasta -out /home/xue/laevis_transcriptome_mar2018/laevis_gg_trancsriptome/blastn_transcriptome_genome/subset_blastout/subset_25_blastout.tsv


```
gmap: 
```
```

### laevis-tropicalis orthologs
### borealis-laevis-tropicalis orthologs


