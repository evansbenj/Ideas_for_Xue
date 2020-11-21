# Raw Data

Raw data were downloaded to the lab computer (~/BenF_GradSchool/Grad_Projects/Borealis_SexChromosome_Family_Transcriptomes/) and to info at (/net/infofile2/2/scratch/benf/Borealis-Family-Transcriptomes-July2017/), which there is a symbolic link to on the cluster in my home directory.

Home directory of data is

*Location:/Borealis-Family-Transcriptomes-July2017/Data/RawReads/*


## Quality assessment

Going to check quality of sequencing with fastqc.

*Location:/Borealis-Family-Transcriptomes-July2017/Analyses/FastQC/*
```bash
screen -S boys -d -m bash -c 'for i in ../../Data/RawReads/B*boy*gz ; do fastqc --threads=4 $i ; done'
screen -S girls -d -m bash -c 'for i in ../../Data/RawReads/B*girl*gz ; do fastqc --threads=4 $i ; done'
screen -S mom -d -m bash -c 'for i in ../../Data/RawReads/B*mom*gz ; do fastqc --threads=4 $i ; done'
screen -S dad -d -m bash -c 'for i in ../../Data/RawReads/B*dad*gz ; do fastqc --threads=4 $i ; done'
```

## Trimmomatic

Testing on one file to make sure the TruSeq is specified correctly.

```bash
screen -S dadLv -d -m java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 BJE3896_dad_liver_R1.fastq.gz BJE3896_dad_liver_R2.fastq.gz BJE3896_dad_liver_R1_paired.fastq.gz BJE3896_dad_liver_R1_unpaired.fastq.gz BJE3896_dad_liver_R2_paired.fastq.gz BJE3896_dad_liver_R2_unpaired.fastq.gz ILLUMINACLIP:/home/benf/bin/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

Will need some semi-intelligent script to pair the files and loop across them all.

*Location:/Borealis-Family-Transcriptomes-July2017/data/RawReads/[mom|dad|girls|boys]/*

```bash
# will operate in a directory, find files and pair them up/created necessary names for Trimmomatic.
/Borealis-Family-Transcriptomes-July2017/Scripts/run-trimmomatic.pl .
```

## Scythe

The fastqc post-Trimmomatic indicates 3' contamination. Run Scythe. (should have done in the reverse order, but the effect should be minimal, just may miss some adapters).

```bash
# Trying in dad dir
for i in *_paired.fastq.gz ; do name=$(grep -o "BJE[0-9]*_[a-z]*_[a-z]*_[A-Z][0-9]" <(echo $i)) ; /home/benf/bin/scythe/scythe -a /home/benf/bin/Trimmomatic-0.36/adapters/trimming-seqs.fa -p 0.1 $i | gzip > $name\_scythe.fastq.gz ; done

# now running on mom/daughters/sons
```

If this works in dad and the quality looks good, then remove the trimmomatic files, as they will not be necessary.

Need to check all with fastqc again to confirm it's all good.

Looked good for dad, some odd kmer stuff hanging around, but pretty minor. Could try trimmomatic again, but likely unnecessary. These odd start kmers are probably related to biases in hexamer priming (which is expected). So I ran on all files.

I then deleted the trimmomatic outputs, keeping just the scythe.
