# Examining the genomic location of DE transcripts 
**Purpose**: When I tried to identify orthologs for borealis DE transcripts, I can only identify orthologs for 10%-25% of the DE transcripts depends on the approachs that I was using. It is odd that so few of them have orthologs in laevis. So I want to examin the DE transcripts themselves and see whether it is due to incorrect assembly.

Path to files that I will use:
```
#path to directory that I will store the output of those analysis
/home/xue/borealis_DE/liver_mvsf/borealis_laevis_tropicalis_orthologs/pre-identification_examination

#path to borealis_De_laevis_genome_blastout
/home/xue/borealis_DE/liver_mvsf/mapping_blastn/liver_mvsf_blastout

```



# Examining the DE_laevis_genome_blastout
The questions that I am interest in looking at is 
  - did transcripts got assembled corrrectly?
  - how many of them were assembled with parts from L&S?
  - Take the one that assembled correctly, extract their genomic location (include all exons, not just hte tophit)
### did transcripts got assembled corrrectly?
To answer the first question, I run the perl script that I developped for my undergraduate thesis with the DE_laevis_genome_blastout
```
cd /home/xue/borealis_DE/liver_mvsf/borealis_laevis_tropicalis_orthologs/pre-identification_examination 

perl ~/script/BlastnOutputAnalysis_4168_Feb20.pl /home/xue/borealis_DE/liver_mvsf/mapping_blastn/liver_mvsf_blastout > /home/xue/borealis_DE/liver_mvsf/borealis_laevis_tropicalis_orthologs/pre-identification_examination/problematic_transcripts.tsv
```
just out of curiosity: count the number of mapped transcripts; the result is 432. The total number of DE transcripts is 434, so it means that 2 DE transcripts was not mapped(why? maybe i can extract their sequence and do a separate blastn mapping). 
```
cd /home/xue/borealis_DE/liver_mvsf/mapping_blastn
cut -f 1 liver_mvsf_blastout | sort -u | wc -l
```
I think I should re-run the blastn to make sure it was done properly when I did it like 7 months ago.
```

```
