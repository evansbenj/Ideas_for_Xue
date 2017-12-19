# Blasting to Trinity-built Borealis Transcriptome
I blasted the below transcript sequence to 6 individual trinity-built Borealis transcriptomes, those were built by BenF.
```
>TRINITY_DN162482_c0_g2_i1
CCCAATATAATCTATATTTTACTGCAAACCATTTCATGCCATTTATATAGTGAATAAAGT
ACCCCTTCTTGGAAAATATAAGCATATTATAAGTTAAAGAGAAGTTCCATGACAATATAA
AAGGACAAGGCCGAAGGCAGAGTGTTTTTATACAGGCCATGGAACTCTGAGGTAACCTCT
AATATCCTCATATTGTACAGCATGGGG
```
The command that I used for constructing the blast database and running blast:
```
makeblastdb -in /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/Individually/BJE3929_boy/BJE3929_boy_trinity.Trinity.fasta -dbtype nucl -out db_BJE3929_boy
makeblastdb -in /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/Individually/BJE4009_girl/BJE4009_girl_trinity.Trinity.fasta -dbtype nucl -out db_BJE4009_girl
makeblastdb -in /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/Individually/BJE4017_boy/BJE4017_boy_trinity.Trinity.fasta -dbtype nucl -out db_BJE4017_boy
makeblastdb -in /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/Individually/BJE4039_boy/BJE4039_boy_trinity.Trinity.fasta -dbtype nucl -out db_BJE4039_boy
makeblastdb -in /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/Individually/BJE4072_girl/BJE4072_girl_trinity.Trinity.fasta -dbtype nucl -out db_BJE4072_girl
makeblastdb -in /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/Individually/BJE4082_girl/BJE4082_girl_trinity.Trinity.fasta -dbtype nucl -out db_BJE4082_girl
makeblastdb -in /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/Individually/Dad/Dad_trinity.Trinity.fasta -dbtype nucl -out db_Dad
makeblastdb -in /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/Individually/Mom/mom_trinity.Trinity.fasta -dbtype nucl -out db_Mom

blastn -task blastn -db /home/xue/borealis_int_seq/transcriptome_blastdb/db_BJE3929_boy -outfmt 6 -evalue 0.00005 -query /home/xue/borealis_int_seq/int_transcript.fasta -out /home/xue/borealis_int_seq/BJE3929_boy
blastn -task blastn -db /home/xue/borealis_int_seq/transcriptome_blastdb/db_BJE4009_girl -outfmt 6 -evalue 0.00005 -query /home/xue/borealis_int_seq/int_transcript.fasta -out /home/xue/borealis_int_seq/BJE4009_girl
blastn -task blastn -db /home/xue/borealis_int_seq/transcriptome_blastdb/db_BJE4039_boy -outfmt 6 -evalue 0.00005 -query /home/xue/borealis_int_seq/int_transcript.fasta -out /home/xue/borealis_int_seq/BJE4017_boy
blastn -task blastn -db /home/xue/borealis_int_seq/transcriptome_blastdb/db_BJE4017_boy -outfmt 6 -evalue 0.00005 -query /home/xue/borealis_int_seq/int_transcript.fasta -out /home/xue/borealis_int_seq/BJE4039_boy
blastn -task blastn -db /home/xue/borealis_int_seq/transcriptome_blastdb/db_BJE4072_girl -outfmt 6 -evalue 0.00005 -query /home/xue/borealis_int_seq/int_transcript.fasta -out /home/xue/borealis_int_seq/BJE4072_girl
blastn -task blastn -db /home/xue/borealis_int_seq/transcriptome_blastdb/db_BJE4082_girl -outfmt 6 -evalue 0.00005 -query /home/xue/borealis_int_seq/int_transcript.fasta -out /home/xue/borealis_int_seq/BJE4082_girl
blastn -task blastn -db /home/xue/borealis_int_seq/transcriptome_blastdb/db_Dad -outfmt 6 -evalue 0.00005 -query /home/xue/borealis_int_seq/int_transcript.fasta -out /home/xue/borealis_int_seq/Dad
blastn -task blastn -db /home/xue/borealis_int_seq/transcriptome_blastdb/db_Mom -outfmt 6 -evalue 0.00005 -query /home/xue/borealis_int_seq/int_transcript.fasta -out /home/xue/borealis_int_seq/Mom
```

Below is the BLAST results in output format6:
```
==> BJE3929_boy <==
TRINITY_DN162482_c0_g2_i1       TRINITY_DN47339_c23_g1_i3       87.879  165     19      1       43      207     346     183     1.61e-51        205
TRINITY_DN162482_c0_g2_i1       TRINITY_DN47339_c23_g1_i3       90.141  71      7       0       43      113     71      1       5.99e-19        96.9
TRINITY_DN162482_c0_g2_i1       TRINITY_DN47339_c23_g1_i3       73.006  163     42      2       46      207     166     327     8.90e-17        89.7
TRINITY_DN162482_c0_g2_i1       TRINITY_DN52217_c0_g4_i2        86.503  163     22      0       43      205     149     311     8.32e-49        196
TRINITY_DN162482_c0_g2_i1       TRINITY_DN52217_c0_g4_i2        72.143  140     37      2       69      207     306     168     8.33e-11        69.8
TRINITY_DN162482_c0_g2_i1       TRINITY_DN58414_c4_g1_i1        85.366  164     23      1       44      207     197     35      1.50e-45        185

==> BJE4009_girl <==
TRINITY_DN162482_c0_g2_i1       TRINITY_DN32568_c2_g4_i5        85.806  155     22      0       40      194     166     12      1.29e-44        181
TRINITY_DN162482_c0_g2_i1       TRINITY_DN32568_c2_g4_i5        72.143  140     37      2       69      207     6       144     5.87e-11        69.8
TRINITY_DN162482_c0_g2_i1       TRINITY_DN32568_c2_g4_i3        85.806  155     22      0       40      194     166     12      1.29e-44        181
TRINITY_DN162482_c0_g2_i1       TRINITY_DN32568_c2_g4_i3        72.143  140     37      2       69      207     6       144     5.87e-11        69.8
TRINITY_DN162482_c0_g2_i1       TRINITY_DN32568_c2_g4_i1        85.806  155     22      0       40      194     166     12      1.29e-44        181
TRINITY_DN162482_c0_g2_i1       TRINITY_DN32568_c2_g4_i1        72.143  140     37      2       69      207     6       144     5.87e-11        69.8

==> BJE4017_boy <==
TRINITY_DN162482_c0_g2_i1       TRINITY_DN46244_c2_g1_i1        85.714  168     23      1       40      207     9       175     1.02e-47        192
TRINITY_DN162482_c0_g2_i1       TRINITY_DN46244_c2_g1_i1        72.393  163     43      2       46      207     192     31      3.82e-15        84.2
TRINITY_DN162482_c0_g2_i1       TRINITY_DN47556_c5_g1_i1        85.976  164     23      0       44      207     494     331     1.02e-47        192
TRINITY_DN162482_c0_g2_i1       TRINITY_DN47556_c5_g1_i1        73.718  156     39      2       53      207     322     476     8.98e-17        89.7
TRINITY_DN162482_c0_g2_i1       TRINITY_DN46244_c2_g1_i17       87.582  153     19      0       42      194     21      173     3.57e-47        190
TRINITY_DN162482_c0_g2_i1       TRINITY_DN46244_c2_g1_i17       82.301  113     20      0       77      189     298     410     2.25e-24        114

==> BJE4039_boy <==
TRINITY_DN162482_c0_g2_i1       TRINITY_DN46244_c2_g1_i1        85.714  168     23      1       40      207     9       175     1.02e-47        192
TRINITY_DN162482_c0_g2_i1       TRINITY_DN46244_c2_g1_i1        72.393  163     43      2       46      207     192     31      3.82e-15        84.2
TRINITY_DN162482_c0_g2_i1       TRINITY_DN47556_c5_g1_i1        85.976  164     23      0       44      207     494     331     1.02e-47        192
TRINITY_DN162482_c0_g2_i1       TRINITY_DN47556_c5_g1_i1        73.718  156     39      2       53      207     322     476     8.98e-17        89.7
TRINITY_DN162482_c0_g2_i1       TRINITY_DN46244_c2_g1_i17       87.582  153     19      0       42      194     21      173     3.57e-47        190
TRINITY_DN162482_c0_g2_i1       TRINITY_DN46244_c2_g1_i17       82.301  113     20      0       77      189     298     410     2.25e-24        114

==> BJE4072_girl <==
TRINITY_DN162482_c0_g2_i1       TRINITY_DN35835_c36_g1_i1       87.654  162     20      0       46      207     265     104     3.75e-51        203
TRINITY_DN162482_c0_g2_i1       TRINITY_DN35835_c36_g1_i1       74.233  163     40      2       46      207     88      249     1.15e-19        98.7
TRINITY_DN162482_c0_g2_i1       TRINITY_DN27653_c1_g1_i2        86.667  165     21      1       43      207     304     467     5.57e-49        196
TRINITY_DN162482_c0_g2_i1       TRINITY_DN27653_c1_g1_i2        70.748  147     40      3       54      199     475     331     2.89e-08        60.8
TRINITY_DN162482_c0_g2_i1       TRINITY_DN32084_c4_g1_i3        85.455  165     24      0       43      207     244     80      2.37e-47        190
TRINITY_DN162482_c0_g2_i1       TRINITY_DN32084_c4_g1_i3        82.895  152     24      2       56      207     492     343     2.08e-35        150

==> BJE4082_girl <==
TRINITY_DN162482_c0_g2_i1       TRINITY_DN33886_c3_g3_i1        89.103  156     17      0       41      196     661     816     1.12e-51        205
TRINITY_DN162482_c0_g2_i1       TRINITY_DN33886_c3_g1_i2        89.103  156     17      0       41      196     1579    1734    1.12e-51        205
TRINITY_DN162482_c0_g2_i1       TRINITY_DN40694_c1_g1_i1        87.117  163     19      1       43      205     106     266     1.66e-49        197
TRINITY_DN162482_c0_g2_i1       TRINITY_DN38339_c3_g1_i5        85.976  164     22      1       44      207     224     386     8.62e-47        188
TRINITY_DN162482_c0_g2_i1       TRINITY_DN38339_c3_g1_i5        69.697  165     45      5       46      207     404     242     1.28e-06        55.4
TRINITY_DN162482_c0_g2_i1       TRINITY_DN38848_c8_g4_i2        87.838  148     18      0       43      190     150     297     3.01e-46        187

```
