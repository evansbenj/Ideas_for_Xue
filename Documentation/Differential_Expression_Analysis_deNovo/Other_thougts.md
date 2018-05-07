# Binning trancripts 
The Perl script that do the binning would: 
  - group overlapping transcripts into bins.  
  - count the number of transcripts in each bin. If 100k transcripts were mapped to chr8L, how many of them overlapped?
  - sum the expression level of transcripts in the same bin -> see if bins on chr8L non-recombinating region have different expression level than other bins (on chr8L non-sex-linked region and other autosomes)
  


# Compare to *X.laevis* chr8L
### *De novo* assembled *X.laevis* transcriptome
The path to *X.laevis* RNA-seq trimmed data and assembled transcriptome: /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Session-SRA. The path to the assembled transcriptome: /home/xue/borealis_DE/Session_Xl_DeNovo_Transcriptome  
```
/home/xue/software/trinityrnaseq-Trinity-v2.5.1/Trinity --seqType fq --left /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Session-SRA/*R1_scythe*  --right /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Session-SRA/*R2_scythe*  --CPU 20 --full_cleanup --max_memory 200G  --output /home/xue/borealis_DE/Session_Xl_DeNovo_Transcriptome/Xl_Transcript_TrinityOut
```

- do the same pipleline with XL RNA-seq data
- compare the expression level of DE genes in XB male, XB female, XL male, and XL female
- do a volcano plot to visualize those expression levels 

# Do trinity assembly for XB transcriptome with different setting
how did others do polyploid transcriptome assembly using Trinity? -> mostly default setting; 
```
Duan, J., Xia, C., Zhao, G., Jia, J. &Kong, X. Optimizing de novo common wheat transcriptome assembly using short-read RNA-Seq data. BMC Genomics 13, 392 (2012).
> ALLPATHSLG error correction, and the paired fragment length was set to 200 bp

Chopra, R. et al. Comparisons of De Novo Transcriptome Assemblers in Diploid and Polyploid Species Using Peanut (Arachis spp.) RNA-Seq Data.
> the default kmer of 25, minimum coverage of 2

Nakasugi, K., Crowhurst, R., Bally, J., Waterhouse, P. &Mittapalli, O. Combining Transcriptome Assemblies from Multiple De Novo Assemblers in the Allo-Tetraploid Plant Nicotiana benthamiana. PLoS One 9, (2014).
> default setting

Gutierrez-Gonzalez, J. J., Tu, Z. J., & Garvin, D. F. (2013). Analysis and annotation of the hexaploid oat seed transcriptome. BMC genomics, 14(1), 471.
> minimum assembled transcript length = 100 nt

Nakasugi, K., Crowhurst, R., Bally, J., Waterhouse, P. &Mittapalli, O. Combining Transcriptome Assemblies from Multiple De Novo Assemblers in the Allo-Tetraploid Plant Nicotiana benthamiana. PLoS One 9, (2014).
> Reads were normalized by the perl script insilico_read_normalization.pl (Trinity package) for Trinity;
> Two settings of k-mer 25 and k-mer 31

He, B. et al. Optimal assembly strategies of transcriptome related to ploidies of eukaryotic organisms. (2011). doi:10.1186/s12864-014-1192-7
> Trinity_release_20131110 (http://trinityrnaseq.github.io/) [16], which used in default parameter, kmer =25, was used in SASP strategy and its Command-line parameters were “--seqType fq --left Reads_1.fq --right Reads_2.fq --CPU 20”

McKain, M. R. et al. A phylogenomic assessment of ancient polyploidy and genome evolution across the Poales. Genome Biol. Evol. 8, evw060 (2016).
> Trinity (Release 2012-06-08; Grabherr et al. 2011 ) with default parameters

Li, H.-Z. et al. Evaluation of Assembly Strategies Using RNA-Seq Data Associated with Grain Development of Wheat (Triticum aestivum L.). PLoS One 8, e83530 (2013).
> Trinity (release 20120608) with minimum contig length  =  100, default parameter

Chen, S., McElroy, J. S., Dane, F. &Peatman, E. Optimizing Transcriptome Assemblies for Leaf and Seedling by Combining Multiple Assemblies from Three De Novo Assemblers. Plant Genome 8, 0 (2015).
> normalized with Trinity’s in silico read normalization (http://Trinityrnaseq.sourceforge.net/), with maximum coverage of 30 and k-mer size of 25
>  The merged assembly was first processed by CD-HIT-EST v4.6.1-2012-08-27 (http://weizhong-lab.ucsd.edu/cd-hit/) with the strictest parameters “-c 1.0-n 10” to remove identical fragments and then subjected to the EvidentialGene tr2aacds pipeline (http://arthropods.eugenes.org/EvidentialGene/about/EvidentialGene_trassembly_pipe.html). 

```
Try #1: increase min contig length to 300 -> stopped and abandoned
- By increasing the filter limit for minimum contig length (--min_contig_length), short fragment that is smaller than the limit would not be output by Trinity. Default is 200 so I increase it to 300. However, this might lead to some information loss.
```
/home/xue/software/trinityrnaseq-Trinity-v2.5.1/Trinity --seqType fq --left /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/*R1*  --right /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/*R2*  --CPU 20 --full_cleanup --max_memory 200G --min_contig_length 300 
```

Try #2: change the kmer size to 31 -> still in progress
- Maximum kmer size for Trinity is 32. One of the above paper set it to 31. 
Accouding to this paper (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5237644/), assembly produces less transcripts with larger kmer size. 
```
/home/xue/software/trinityrnaseq-Trinity-v2.5.1/Trinity --seqType fq --left /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/*R1*  --right /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/*R2*  --CPU 20 --full_cleanup --max_memory 200G --min_kmer_cov 31 --output /home/xue/xborealis_alltogether/xb_kmer31_trinity
```

Try #3: with RNA-seq from liver only -> stopped and abandoned 
- Since sex-specific genes could have different expression level in liver and gonad. When doing expression analysis, 
```
/home/xue/software/trinityrnaseq-Trinity-v2.5.1/Trinity --seqType fq --left /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/*liver_R1*  --right /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/*liver_R2*  --CPU 20 --full_cleanup --max_memory 200G  --output /home/xue/xborealis_alltogether/xb_liver_trinity
```
Error message:
```
Error: Encountered internal Bowtie 2 exception (#1)
Command: bowtie2-build --wrapper basic-0 --threads 20 -o 3 /home/xue/xborealis_alltogether/xb_liver_trinity/chrysalis/inchworm.K25.L25.DS.fa.min100 /home/xue/xborealis_alltogether/xb_liver_trinity/chrysalis/inchworm.K25.L25.DS.fa.min100
Error, cmd: bowtie2-build --threads 20 -o 3 /home/xue/xborealis_alltogether/xb_liver_trinity/chrysalis/inchworm.K25.L25.DS.fa.min100 /home/xue/xborealis_alltogether/xb_liver_trinity/chrysalis/inchworm.K25.L25.DS.fa.min100 1>/dev/null 2>tmp.39623.stderr died with ret 256 at /home/xue/software/trinityrnaseq-Trinity-v2.5.1/PerlLib/Pipeliner.pm line 166.
        Pipeliner::run(Pipeliner=HASH(0x1f440d0)) called at /home/xue/software/trinityrnaseq-Trinity-v2.5.1/Trinity line 1827
        main::run_chrysalis("/home/xue/xborealis_alltogether/xb_liver_trinity/inchworm.K25"..., "/home/xue/xborealis_alltogether/xb_liver_trinity/both.fa", 200, 500, undef, "/home/xue/xborealis_alltogether/xb_liver_trinity/both.fa", "/home/xue/xborealis_alltogether/xb_liver_trinity/both.fa") called at /home/xue/software/trinityrnaseq-Trinity-v2.5.1/Trinity line 1667
        main::run_Trinity() called at /home/xue/software/trinityrnaseq-Trinity-v2.5.1/Trinity line 1320
        eval {...} called at /home/xue/software/trinityrnaseq-Trinity-v2.5.1/Trinity line 1319

Trinity run failed. Must investigate error above.
```
# things that we could do but not urgent
- repeat DE analysis with different softwware package (RESM, DESeq2 etc)
- repeat the pipeline with DE comparison pair (all male smaples vs all female samples, male gonad vs female gonad, liver vs gonad)

# Path to session laevis RNAseq read data
```
/home/benf/Borealis-Family-Transcriptomes-July2017/Data/Laevis-Session-SRA
```

