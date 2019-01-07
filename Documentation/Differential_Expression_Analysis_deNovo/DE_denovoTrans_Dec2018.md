# Repeat DE for borealis transcriptome Dec2018
## transcriptome building
Graham has the newest version of Trinity, hence, I rebuild the transcriptome with the following bash script. 
```

```

## Mapping borealis transcriptome to laevis genome
```
time gmap -D /home/songxy/scratch/tropicalis_transcriptome/tropicalis_genome/db_gmap_tropicalis_v91 -d db_gmap_tropicalis_v91 -A -B 5 -t 15 -f samse /home/songxy/scratch/tropicalis_transcriptome/transcriptome_building/tropicalis_transcriptome_trinityOut.Trinity.SuperTrans.fasta | samtools view -S -b > /home/songxy/scratch/tropicalis_transcriptome/mapping_transcriptome_to_genome/tropicalis_denovoT_tropicalisv91_genome_gmap.bam
```
