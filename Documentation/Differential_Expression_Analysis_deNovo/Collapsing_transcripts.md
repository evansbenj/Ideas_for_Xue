# Callapsing transcripts
Collapsing the 1.6 million transcripts to a number that is more biologically realistic (i.e. remove all of the false transcripts due to sequence error, polymorphism, etc..).

The callapsing of transcripts was done using a perl script. The script can be found in the script folder on info (path: /home/xue/script/identify_genomic_location_for_transcripts.pl) . Here I will illustrate the basic sctructure of the script. 

 - This perl script will take in the mapping file, in which the transcriptome was mapped to a genome and will provide the genomic location of mapping.
 - With the aid of the laevis gff file, which contain genomic coordinates of all annotated genes, the script can figure out which gene the transcripts mapped to based on genomic coordinate.
 - Afterward, the script will read in the expression value of each transcript (output of kallisto) and sum the expression value of transcripts that mapped to the same gene.
 - Lastly, it will output the callapsed expression value for each gene.


#### generate gene only gff file 
extract gene information from gff file
```
awk '$3 == "gene"{print}' XENLA_9.2_Xenbase.gff3 > /home/xue/genome_data/laevis_genome/gff3_xl9_2/XENLA_gene.gff3
```

##### generate the alignment file
mapping transcriptome to genome
```
time gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92/ -d laevis92_gmap -A -B 5 -t 8 -f samse --cross-species /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/subset_1.fasta | samtools view -S -b > /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/borealis_denovoT_laevisV92_genome_gmap_subset1.bam
```
filter the mapping read remove the unmapped reads, which is indicated by flag 0x04. Then I piled the filtered output to bedtools, which extracts alignment coordinates for transcripts based on their CIGAR strings;
```
samtools view -F 0x04 -b borealis_denovoT_laevisV92_genome_gmap.bam | bedtools bamtobed -i > borealis_denovoT_laevisV92_genome_gmap_bedfile.bed
```
filter the bedfile and keep only the hit with the highest mapping quality score; if there are two hit with the same mapping quality score, keep the first one only;
```
perl ~/script/orthologs_identification_filterBedfile.pl borealis_denovoT_laevisV92_genome_gmap_bedfile.bed > borealis_de_laevisV92_genome_gmap_bedfile_filtered.tsv
```
##### expression matrix 
```
#Indexing the transcriptomes:
kallisto index -i /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta.kallisto_idx /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta
# Transcript abundance quantification:
kallisto quant -i /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/trinity_out_dir.Trinity.fasta.kallisto_idx  -o female_rep7 <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4082_girl_liver_R1_scythe.fastq.gz) <(gunzip -c /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE4082_girl_liver_R2_scythe.fastq.gz)
# compile the read count for each samples into a matrix; used a script included in the Trinity RNAseq analysis package
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/util/abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix borealis_liver  --name_sample_by_basedir female_rep1/abundance.tsv female_rep2/abundance.tsv female_rep3/abundance.tsv female_rep4/abundance.tsv male_rep1/abundance.tsv male_rep2/abundance.tsv male_rep3/abundance.tsv male_rep4/abundance.tsv
```

#### Path to all my input files
```
#laevis gff
#####/home/xue/genome_data/laevis_genome/gff3_xl9_2/XENLA_gene.gff3
#alignment files: borealis, laevis, tropicals
### /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/borealis_denovoT_laevisV92_genome_gmap_filtered.bed
/home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_laevisGenomeApproach/laevis_denovoT_laevisV92_genome_gmap_bedfile_filtered.bed
/home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_laevisGenomeApproach/tropicalis_denovoT_laevisV92_genome_gmap_bedfile_filtered.bed

#expression matrix: borealis, laevis, tropicalis
/home/xue/borealis_DE/de_sex_liver/kallisto_out/borealis_liver.counts.matrix
/home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/kallisto_out/laevis_denovo.counts.matrix
/home/xue/tropicalis_transcriptome/tropicalis_denovo_transcriptome/tropicalis_kallisto_denovo/tropicalis_denovo.counts.matrix
```

Path to all the output files
```
/home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/sum_expression/borealis_expression_matrix_per_gene.tsv

/home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/sum_expression/laevis_expression_matrix_per_gene.tsv

/home/xue/tropicalis_transcriptome/tropicalis_denovo_transcriptome/sum_expression/tropicalis_expression_matrix_per_gene.tsv
```

My usage of the perl script
```
time perl ~/script/identify_genomic_location_for_transcripts.pl

time perl ~/script/identify_genomic_location_for_transcripts.pl --align /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_laevisGenomeApproach/laevis_denovoT_laevisV92_genome_gmap_bedfile_filtered.bed --matrix /home/xue/laevis_transcriptome_mar2018/laevis_denovo_transcriptome/kallisto_out/laevis_denovo.counts.matrix > laevis_expression_matrix_per_gene.tsv

time perl ~/script/identify_genomic_location_for_transcripts.pl --align /home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_laevisGenomeApproach/tropicalis_denovoT_laevisV92_genome_gmap_bedfile_filtered.bed --matrix /home/xue/tropicalis_transcriptome/tropicalis_denovo_transcriptome/tropicalis_kallisto_denovo/tropicalis_denovo.counts.matrix > tropicalis_expression_matrix_per_gene.tsv
```

An example of the output format
```
        female_rep1     female_rep2     female_rep3     female_rep4     male_rep1     male_rep2     male_rep3     male_rep4

Scaffold1000:541-2366   30      21.99997        27.00004        32.00002
        69      19.99997
Scaffold100:209407-303984       0       0       1.56265 0       0       2
Scaffold100:35531-36473 1       0       0       0       0       1
Scaffold10166:48-329    0       0       2       0       0       1
Scaffold101:16971-37533 51.55615        54.30357        47.3974 103.02645
       56.63398        110.96152
Scaffold101:184491-195741       6       2       2       2       2       4
Scaffold101:60857-77058 1       0       1       1       0       0
Scaffold1024:84-3527    50.910796       15.27986        78.06537        25.40181        32.01005        114.97102
Scaffold102:49146-137487        70.3708 113.54159       80.1207 123.97537
       99.96169        190.05188
Scaffold103:12613-17294 6.79133 8.1832  4       25.7545 2       11.0367
Scaffold103:126312-133620       2       0       6.48314 2.13222 7       7.60564
Scaffold103:141280-147297       736.9202        638.664019      621.780150289   1143.13413694   1019.712        1103.796
Scaffold103:154382-166533       892.44722       822.81503       729.29505
       1310.92532      1181.29372      1486.3976
Scaffold103:182117-199709       19      7       12      22.1954 33.2024 14
Scaffold103:213204-217315       1       0       2       4       3       1
Scaffold103:223450-227723       0       0       0       1       0       1
Scaffold103:245052-268283       205.6287        135.01758       147.70356
       214.66176       233.65407       231.27011
Scaffold103:26413-43870 1250.89195      1099.32262      1328.82234      1921.47163      2110.81894      1856.720659
Scaffold103:279114-285495       6       0       7       3       7       10.1566
Scaffold103:292162-300280       9.460585        0       5       5.31976 4.25588 10.28439
Scaffold103:49574-80934 11.81191        11.75369        16.48935        15.50125        12.47701        20.87713
Scaffold104:120562-221333       25.00002        21.8177 100.07563       73.71639        74.99995        58.55051
Scaffold104:283726-298386       616.74527       405.816 569.33025       984.4[xue@info115 ~]$ less borealis_transcriptome/borealis_denovo_transcriptome_august2017/sum_expression/borealis_expression_matrix_per_gene.tsv
         female_rep1     female_rep2     female_rep3     female_rep4     male_rep1
MT:10711-11491  973972  299235  179192  227215  276145  327805  229833  375904
MT:16249-17388  561066  236644  139493  151080  129554  333908  206843  252375
Scaffold1000:541-2366   0       0       1       0       0       0       0
       0
Scaffold100:110725-181202       0       11      9       6       4       10
      21      1
Scaffold100:186956-193508       7       1       6       2       2.41232 1
       0       5.76162
Scaffold100:209407-303984       43      105.5498        58      73      98
      33.20251        136     40.0171
Scaffold100:52396-54980 0       0       0       0       0       1       0
       0
Scaffold100:71321-100665        0       2       3       0       0       0
       0       0
Scaffold100:85642-93754 0       2       2       0       0       0       0
       0
Scaffold101:16971-37533 0       4.59942 3.04901 6.05744 3.02576 3.05254 5
       3.13822
Scaffold101:184491-195741       1       0       0       1       0       0
       0       0
Scaffold101:60857-77058 0       1       0.293095        1       1       0
       1       1.895972
Scaffold1022:207-1128   128.56291       73.619661859    114.72426       72.054  76.67174        81.44043        17.9909 54.6137
Scaffold102:289198-305512       36.3282 3.47788 1.2831  1.481656        3.85323 17.9857 5.86822 2
Scaffold102:49146-137487        66      83.48474        44.84453        49.12221        25.38133        45.82277        42.79529        51.81923
Scaffold102:90556-91596 0       1       0       0       2       0       0
       0
         female_rep1     female_rep2     female_rep3     female_rep4     maleMT:10711-11491  973972  299235  179192  227215  276145  327805  229833  37590MT:16249-17388  561066  236644  139493  151080  129554  333908  206843  25237Scaffold1000:541-2366   0       0       1       0       0       0       0
Scaffold100:110725-181202       0       11      9       6       4       10
Scaffold100:186956-193508       7       1       6       2       2.41232 1    Scaffold100:209407-303984       43      105.5498        58      73      98
         female_rep1     female_rep2     female_rep3     female_rep4     maleMT:10711-11491  973972  299235  179192  227215  276145  327805  229833  37590MT:16249-17388  561066  236644  139493  151080  129554  333908  206843  25237Scaffold1000:541-2366   0       0       1       0       0       0       0
Scaffold100:110725-181202       0       11      9       6       4       10
Scaffold100:186956-193508       7       1       6       2       2.41232 1
Scaffold100:209407-303984       43      105.5498        58      73      98   Scaffold100:52396-54980 0       0       0       0       0       1       0    Scaffold100:71321-100665        0       2       3       0       0       0    Scaffold100:85642-93754 0       2       2       0       0       0       0
Scaffold101:16971-37533 0       4.59942 3.04901 6.05744 3.02576 3.05254 5
Scaffold101:184491-195741       1       0       0       1       0       0    Scaffold101:60857-77058 0       1       0.293095        1       1       0
Scaffold1022:207-1128   128.56291       73.619661859    114.72426       72.05Scaffold102:289198-305512       36.3282 3.47788 1.2831  1.481656        3.853Scaffold102:49146-137487        66      83.48474        44.84453        49.12Scaffold102:90556-91596 0       1       0       0       2       0       0
Scaffold1032:667-1734   0       2       0       0       0       0       0
Scaffold103:12613-17294 25.29398        31.339  26      25.87521        27.32Scaffold103:126312-133620       3       1       1       2       1       5    Scaffold103:141280-147297       49.40628        160.08408       75.512197
Scaffold103:154382-166533       333.7233        1097.80433      963.22563    Scaffold103:182117-199709       23.1372 3       1.999998036     3.62443 2    Scaffold103:213204-217315       20.959912       13.575723       3.509758     Scaffold103:223450-227723       0       0       0       0       0       0
Scaffold103:245052-268283       140.28373       158.8871        193.0601
Scaffold103:26413-43870 286.576 356.927 193.405 210.188 0       221.695 296.2Scaffold103:279114-285495       1       0       0       0       0       0
Scaffold103:292162-300280       19.3758 14.0796 4.3328  0       22.5663 11.99Scaffold103:3229-4058   0       31.4519 9.92041 13.165  17.2759 8.00091 38.96Scaffold103:49574-80934 276.81823       153.16803       92.760987       211.4Scaffold1045:1830-2841  5       1       0       4       3       1       2
Scaffold104:109217-111925       377.9995        6       5       2       16   Scaffold104:120562-221333       4       7       5.01495 9       5       2
Scaffold104:208710-216564       0       0       0       0       1       0
Scaffold104:283726-298386       318.557748036887        446.700944      327.9Scaffold104:68492-71626 0       1       0       1       0       2       0    Scaffold104:74556-83569 53.3479 13.7573 8.52446 16.80694        26.9081592
Scaffold104:83543-100129        45.99998        50.18748        37.96139
:
```


