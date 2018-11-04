# Building Transcriptomes

Trinity manual suggests doing a single build across all samples for DE analysis, and I guess that makes sense. Afterwards, DE analysis involves mapping reads back to a transcriptome to get level counts, thus, either we map to a chosen transcriptome (yuck), or to a single build.

Trinity manual says 1Gb ~ 1Mil reads, the trimmed files have about 53Mil per sample (for each direction...only kept paired reads). So that is 50Gb per sample (or is it 100Gb? they don't say if this is per direction), with 15 samples is at least 750 Gb (and they say 0.5 - 1 h per 1Mil...which is potentially 750 hours, i.e., a month! -- better not be). The builds did not take this much memory or time.

Would individual transcriptomes be useful for sex specific sequences? I guess mapping to a single build would inform us of this as well (e.g., all females map, males do not).


#### Tips
Trinity manual suggest --normalize_reads option for large datasets to improve run time, does this apply to our data?

From manual: "you can reduce the total RAM requirements by running Trinity with parameter '--min_kmer_cov 2'. Although the assembly should still be of high quality and require less RAM, lowly expressed transcripts may be more highly fragmented in the assembly."

* not ideal, as we are looking for things highly expressed in one sex, low in the other. Then again, lacking detection in one sex across all samples is indicative of low expression (but can't tell that from just missing the transcripts).


## Building

Ok, merging all of the files, will transfer to sharcnet (just for memory), and try and do a single build. Will need to confirm the Trinity version Ben has, or just local install my own.

Running on goblin.

```bash
sqsub -q threaded -n 10 --mpp=500g -r 120h -o borealis-allTogether-goblin500GB-120h-clusterOut.txt ./Trinity-Build/Trinity --seqType fq --left borealis_allSamples_left.fastq.gz --right borealis_allSamples_right.fastq.gz --CPU 10 --full_cleanup --max_memory 500G
```

Trying on Info too

*Location:/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/All-together/*

```bash
/home/benf/bin/Trinity-v2.4.0/Trinity --seqType fq --left ../../Trimmed/borealis_allSamples_left.fastq.gz --right ../../Trimmed/borealis_allSamples_right.fastq.gz --CPU 20 --full_cleanup --max_memory 200G

# mapping to the laevis genome, we'll need that eventually. Or could wait and just map transcripts of interest....hmm, this is something to think about b/c exons...GMAP, BLAT...maybe just use mapped reads?
```
So the info build seems to have finished. I'll let the goblin one continue for now, but it's been qued for days. Killed. Need to remove files (password wasn't working).


## Building Individually

I need to do individual builds, maybe, to try and resolve W and Z copies, though they should be present in the complete build too. Minimally, I should do mom and dad builds, then maybe I could do all daughters and all sons as two more builds (again, offspring should get one allele from each and perhaps with more repetition there may be more confidence in each allele construction).


*Location:/Borealis-Family-Transcriptomes-July2017/Data/Trinity-Build-Info/Individually/*
```bash
# mom
/home/benf/bin/Trinity-v2.4.0/Trinity --seqType fq --left ../../Trimmed/BJE3897_mom_liver_R1_scythe.fastq.gz --right ../../Trimmed/BJE3897_mom_liver_R2_scythe.fastq.gz --CPU 15 --full_cleanup --max_memory 200G --output mom_trinity

# dad
cat ../../Trimmed/BJE3896_dad*R1_scythe.fastq.gz > ../../Trimmed/dad_all_left.fastq.gz ; cat ../../Trimmed/BJE3896_dad*R2_scythe.fastq.gz > ../../Trimmed/dad_all_right.fastq.gz

/home/benf/bin/Trinity-v2.4.0/Trinity --seqType fq --left ../../../Trimmed/dad_all_left.fastq.gz --right ../../../Trimmed/dad_all_right.fastq.gz --CPU 15 --full_cleanup --max_memory 200G

# right, do need the offspring individually as daughters could get either Z. Building them altogether could confuse the two alleles. IF these were population samples, could actually build trees, females Z should group with all of the male sequences, and W should be a clade of just females.

# BJE4009_girl
cat ../../../Trimmed/BJE4009_girl*R1_scythe.fastq.gz > ../../../Trimmed/BJE4009_girl_all_left.fastq.gz ; cat ../../../Trimmed/BJE4009_girl*R2_scythe.fastq.gz > ../../../Trimmed/BJE4009_girl_all_right.fastq.gz ; /home/benf/bin/Trinity-v2.4.0/Trinity --seqType fq --left ../../../Trimmed/BJE4009_girl_all_left.fastq.gz --right ../../../Trimmed/BJE4009_girl_all_right.fastq.gz --CPU 15 --full_cleanup --max_memory 200G --output BJE4009_girl_trinity

# BJE4017_boy
cat ../../../Trimmed/BJE4017_boy*R1_scythe.fastq.gz > ../../../Trimmed/BJE4017_boy_all_left.fastq.gz ; cat ../../../Trimmed/BJE4017_boy*R2_scythe.fastq.gz > ../../../Trimmed/BJE4017_boy_all_right.fastq.gz ; /home/benf/bin/Trinity-v2.4.0/Trinity --seqType fq --left ../../../Trimmed/BJE4017_boy_all_left.fastq.gz --right ../../../Trimmed/BJE4017_boy_all_right.fastq.gz --CPU 15 --full_cleanup --max_memory 200G --output BJE4017_boy_trinity

# BJE4039_boy
cat ../../../Trimmed/BJE4039_boy*R1_scythe.fastq.gz > ../../../Trimmed/BJE4039_boy_all_left.fastq.gz ; cat ../../../Trimmed/BJE4039_boy*R2_scythe.fastq.gz > ../../../Trimmed/BJE4039_boy_all_right.fastq.gz ; /home/benf/bin/Trinity-v2.4.0/Trinity --seqType fq --left ../../../Trimmed/BJE4039_boy_all_left.fastq.gz --right ../../../Trimmed/BJE4039_boy_all_right.fastq.gz --CPU 15 --full_cleanup --max_memory 200G --output BJE4039_boy_trinity

# BJE4072_girl
cat ../../../Trimmed/BJE4072_girl*R1_scythe.fastq.gz > ../../../Trimmed/BJE4072_girl_all_left.fastq.gz ; cat ../../../Trimmed/BJE4072_girl*R2_scythe.fastq.gz > ../../../Trimmed/BJE4072_girl_all_right.fastq.gz ; /home/benf/bin/Trinity-v2.4.0/Trinity --seqType fq --left ../../../Trimmed/BJE4072_girl_all_left.fastq.gz --right ../../../Trimmed/BJE4072_girl_all_right.fastq.gz --CPU 15 --full_cleanup --max_memory 200G --output BJE4072_girl_trinity

# BJE3929_boy
cat ../../../Trimmed/BJE3929_boy*R1_scythe.fastq.gz > ../../../Trimmed/BJE3929_boy_all_left.fastq.gz ; cat ../../../Trimmed/BJE3929_boy*R2_scythe.fastq.gz > ../../../Trimmed/BJE3929_boy_all_right.fastq.gz ; /home/benf/bin/Trinity-v2.4.0/Trinity --seqType fq --left ../../../Trimmed/BJE3929_boy_all_left.fastq.gz --right ../../../Trimmed/BJE3929_boy_all_right.fastq.gz --CPU 15 --full_cleanup --max_memory 200G --output BJE3929_boy_trinity


# BJE4082_girl
cat ../../../Trimmed/BJE4082_girl*R1_scythe.fastq.gz > ../../../Trimmed/BJE4082_girl_all_left.fastq.gz ; cat ../../../Trimmed/BJE4082_girl*R2_scythe.fastq.gz > ../../../Trimmed/BJE4082_girl_all_right.fastq.gz ; /home/benf/bin/Trinity-v2.4.0/Trinity --seqType fq --left ../../../Trimmed/BJE4082_girl_all_left.fastq.gz --right ../../../Trimmed/BJE4082_girl_all_right.fastq.gz --CPU 15 --full_cleanup --max_memory 150G --output BJE4082_girl_trinity

```

# Re-Building borealise transcriptome with different parameters
Ian and Ben's suggestion: try to trim raw read with a different trimmomatic paramet to see if we can cut down more low quality sequence. and then build the transcriptome again to see if we can reduce the amount of isoforms
### Trimmomatic
A perl script that was written by BenF was modified to run Trimmomatic on each inidividual rnaseq read files
```bash
perl run_trimmomatic.pl 
```
The Trimmomatic command in the script is:
```bash
java -jar /home/xue/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 $r1_name $r2_name $r1_par_out $r1_unpar_out $r2_par_out $r2_unpar_out ILLUMINACLIP:/home/xue/software/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MAXINFO:30:0.7 MINLEN:36
```
### Trinity de novo transcriptome assembly with version 2.4.0
I run trinity v2.4.0 with an additional parameter `--min_kmer_cov 2` (suggested by Ian), which should reduce the amount of low coverage transcript. 
```bash
time /home/xue/software/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq  --left /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_oct2018/Trimmed/borealis_R1_paired.fastq.gz --right /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_oct2018/Trimmed/borealis_R2_paired.fastq.gz --CPU 20 --full_cleanup --max_memory 200G --min_kmer_cov 2 --output /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_oct2018/; echo "trinity is done at info114 in screen assembly" | mail songxy2@mcmaster.ca
```
The amount of transcripts that we end up having after running the above assembly is 1.5m. The initial run in Sept 2017 had 1.6m. This Oct2018 run have 0.1m less transcripts but still not very good. 

### Trinity de novo transcriptome assembly with version 2.7.0-PRERELEASE
I would like to try the newer version of trinity. I tried to build the newest version v2.8 but it needs the most updated version of cmake 2.9 (current version on info is v2.8.12.2). To install cmake 2.9, it needs g++ v4.7 and current version is v4.4.7. I can't update g++ compilor and hence give up on installing Trinity v2.8. I go down the Trinty version list and install Trinity v2.7.0. Trinity v2.7.0 need a newer version of jellyfish and bowtie2 than what we have on info. Hence, I locally install them and export their path. 
```bash
# installing Trinity v2.7.0
# installation instruction: https://github.com/trinityrnaseq/trinityrnaseq/wiki/Installing-Trinity
make 
make plugins

# to test Trinity
cd sample_data/test_Trinity_Assembly/
./runMe.sh

# installing jellyfish 2.2.5; I tried to install the newest jellyfish v2.2.10 but it needs g++ v4.7 as well
# manual: https://github.com/gmarcais/Jellyfish
./configure --prefix=/home/xue/software/jellyfish-2.2.5/
make -j 4

# installing bowtie2 ()
# manual: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2
make static-libs
make STATIC_BUILD=1
make

# export path
export PATH="/home/xue/software/jellyfish-2.2.5/bin:$PATH"
export PATH="/home/xue/software/bowtie2-2.3.4.3:$PATH"
```
In this assembly run, I included `--no_salmon` because the current version of Salmon on info is 0.8, and trinity2.7 need min Salmon v0.9.1. In order to install salmon v0.9, it needs g++ v4.7 and current version is v4.4.7. I can't update g++ compilor and hence give up on installing salmon. Salmon was used for quick expression estimates and filtering of likely artifacts. By running the assembly with `--no_salmon`, we would lose the opportunity to ultilize the new filtering function of trinity, but I guess it would be fine.  

```bash
time /home/xue/software/trinityrnaseq-Trinity-v2.7.0-PRERELEASE/Trinity --seqType fq --samples_file /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_oct2018/Trimmed/borealis_rnaseq_trimmed_samples_file.tsv  --CPU 25 --inchworm_cpu 15 --full_cleanup --max_memory 200G --min_kmer_cov 2 --no_salmon --output /home/xue/borealis_transcriptome/borealis_denovo_transcriptome_oct2018/borealis_denovo_nov2018trinityOut; echo "trinity is done at info113 in screen assembly" | mail songxy2@mcmaster.ca
```
