# Building Transcriptomes

Trinity manual suggests doing a single build across all samples for DE analysis, and I guess that makes sense. Afterwards, DE analysis involves mapping reads back to a transcriptome to get level counts, thus, either we map to a chosen transcriptome (yuck), or to a single build.

Trinity manual says 1Gb ~ 1Mil reads, the trimmed files have about 53Mil per sample (for each direction...only kept paired reads). So that is 50Gb per sample (or is it 100Gb? they don't say if this is per direction), with 15 samples is at least 750 Gb (and they say 0.5 - 1 h per 1Mil...which is potentially 750 hours, i.e., a month! -- better not be).

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
So the info build seems to have finished. I'll let the goblin one continue for now, but it's been qued for days.


## Building Individually

I need to do individual builds, maybe, to try and resolve W and Z copies, though they should be present in the complete build too. Minimally, I should do mom and dad builds, then maybe I could do all daughters and all sons as two more builds (again, offspring should get one allele from each and perhaps with more repetition there may be more confidence in each allele construction).


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

```
