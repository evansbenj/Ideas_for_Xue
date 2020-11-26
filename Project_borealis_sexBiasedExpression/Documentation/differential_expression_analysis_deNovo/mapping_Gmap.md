
# Identify transcript's genomic location

## Mapping borealis transcriptome to laevis genome
I did it on info like this:
```bash
#the indexing of X. laevis genome was previously done
time gmap -D /home/xue/genome_data/laevis_genome/db_gmap_xl92 -d laevis92_gmap -A -B 5 -t 25 -f samse --cross-species /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_dec2018/transcriptome/borealis_transcriptome_trinityOut.fasta | samtools view -S -b > /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_dec2018/analysis/transcriptome/mapping_xb_denovoTrans_xl_genomev92_gmap/borealis_denovoT_laevisv92_genome_gmap.bam
```
Meanwhile, I also tried it in graham just to see which one finish faster and I would know which I should use if I want a faster result.
```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=10:00:00
#SBATCH --mem=15G
#SBATCH --job-name=borealis_gmap_indexing
#SBATCH --account=def-ben

module load nixpkgs/16.09  gcc/7.3.0
module load gmap-gsnap/2018-07-04
module load samtools/1.9

time gmap -D /home/songxy/projects/def-ben/songxy/genome/laevis_genome/db_gmap_xl92 -d laevis92_gmap -A -B 5 -t 20 -f samse --cross-species /home/songxy/projects/def-ben/songxy/borealis_transcriptome/data/borealis_denovo_transcriptome_dec2018_trinityOut/borealis_transcriptome_trinityOut.Trinity.fasta | samtools view -S -b > /home/songxy/projects/def-ben/songxy/borealis_transcriptome/analysis/mapping_xb_denovoTrans_xl_genomev92_gmap/borealis_denovoT_laevisv92_genome_gmap.bam
```
Gmap on graham ran for 14hrs and the one on infor was still running. Terminated the mappping on infoand transfer the mapping file from graham to info for downstream processing since I found it was easier to do stuffs on info.  

## Coverting to bedfile
Following the pipeline that I did in `Collapsing_transcripts.md`.

I removed the unmapped reads, which is indicated by flag 0x04. Then I piled the filtered output to bedtools, which extracts alignment coordinates for transcripts based on their CIGAR strings;
```
samtools view -F 0x04 -b borealis_denovoT_laevisv92_genome_gmap.bam | bedtools bamtobed -i > borealis_denovoT_laevisV92_genome_gmap_bedfile.bed
```
With a perl script, I filterred the bedfile and keep only the hit with the highest mapping quality score; if there are two hit with the same mapping quality score, keep the first one only; **This might cause problem in the future**. 
  - one of the reason that the transcript mapped to multiple location with the same mapping quality is that the transcript was not assembled correctly. This transcript should be drop from the data set
  - another reason could be that it was mapping to the homologous chromosome. This can be fix by looking at the supertranscriptome and see which chromosome the corresponding supertranscript mapped to. -> **look at it in R** -> this can be used as a extra step to confirm genomic location in the case of doubt.
  
```
perl ~/script/orthologs_identification_filterBedfile.pl borealis_denovoT_laevisV92_genome_gmap_bedfile.bed > borealis_de_laevisV92_genome_gmap_bedfile_filtered.tsv
```
The path to bedfile is 
```
/home/xue/borealis_transcriptome/borealis_denovo_transcriptome_dec2018/analysis/transcriptome/mapping_xb_denovoTrans_xl_genomev92_gmap/
```

