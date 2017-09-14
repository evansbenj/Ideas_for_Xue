# Sorting W and Z copies

Want to figure out the W and Z specific sequences. Will be many transcripts/alleles per gene in the sex linked region. Building trees may help as female W should form a clade, and female Z should form another clade with all male sequences. Could then write a script looking for all female clades, sort of thing.

## Get laevis genes from sex linked region

Need to parse apart the genome using the gff file. Might go with the tool [gffread](http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread_ex), but test what it does with overlapping sequences, which exist in gff files. Could also look at the bedtools getfasta tool.

Ok, confirmed that gffread just grabs the first instance of a gene in a gff. That's not good, as there could be a longer transcript. I could write a script to reduce the gff to just the longest transcripts...then run this. That's probably best, anything I write will be substantially slower.

laevis super scaffold ref here: */home/benf/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/*

Ok, just doing a raw run to see what happens.
```bash
# grab genes in the sex linked region only
grep "chr8L" XL_9.1_v1.8.3.2.primaryTranscripts.gff3 | awk '$5 < 57000000' > chr8L_sub57mil.gff

~/bin/gff-read/gffread/gffread chr8L_sub57mil.gff -g /home/benf/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/Xla.v91.superscaffold.fa -w chr8L_raw_gffread.fa
```

Looks like I may not need to worry about longest as there is the same number of sequences output as there are genes in the gff file (i.e., no genes missing).

```bash
grep -c "gene" chr8L_sub57mil.gff
grep -c ">" chr8L_raw_gffread.fa
```

Would be good to confirm that a few genes that should be there are and that they have the full sequence they should.

For now I have just grabbed < 57000000, which is not the precise boundary of sex linkage according to the cross. So, I'll have to keep inheritance in mind and what not (at least confirm that these three individiuals of each sex have the same seq, which trees should resolve).


Make a blastdb of the chr8L sequences and confirm genes are there. Should probably compare to genome sequence exons+primes.

```bash
makeblastdb -dbtype nucl -in ../chr8L_raw_gffread.fa -out chr8L_raw_gffread

blastn -task blastn -db blast_dbs/chr8L_raw_gffread -query rand-confirmation-gene2.fa -outfmt 6 -evalue 0.0005
```

Checked sox3 and 2 random genes. They were there in the entirety. Looks like the gffread program worked.

## Grabbing borealis transcripts for each gene

Now I need to associate the borealis transcripts with *X. laevis* genes from this region. Will also eventually want *S. tropicalis* and the chr8S sequence. What is my game plan here? Somehow I need to tease apart the 2 alleles, Z and W, from the females. I obviously can't inspect each gene by eye (though there are only ~1000). It would also be good to have something quantitative.

Perhaps an approach similar to the homeolog ID from our 2016 paper with tree building and parsing. The expectation is that females should form 2 clades, one for each allele, so long as there is no recombination between the Z and W. One clade, the Z obviously, should group with the males, and the other should be on it's own. This would work for population samples or cross data. I could include the parents in this and mom should then have a sequence in each clade (presuming that Trinity assembled the alleles correctly). This approach is attractive as there could be many transcripts per gene, and all we are looking for here are two clades (one female specific, the other mixed sex).

An alternative approach would be to collect the offspring transcripts for a particular gene, and then inspect the maternal transcripts for that gene, and track site inheritance. Basically, what I was doing by eye for the sox3 and other genes. Could be interesting coding. With many transcripts per gene, this could get messy to know which ones to compare (between mom and offspring). I'm thinking trees are easier. I could get a clade of W sequences, and pick the longest for each offspring, or something.

Let's test this on a few known genes: sox3, nr5a1 (SF-1), and FMR1.

1. need to make blast dbs of each individual transcriptome
*Location:/home/benf/Borealis-Family-Transcriptomes-July2017/Analyses/ZW_Gene_Seqs/blast_dbs/*
```bash
for i in ../../../Data/Trinity-Build-Info/Individually/*/*fasta ; do echo $i ; dirtname=$(grep -o "BJE[0-9]*_[a-z]*\/" <(echo $i)) ; name=$(grep -o "BJE[0-9]*_[a-z]*" <(echo $dirtname)) ; screen -m -d bash -c "makeblastdb -in $i -dbtype nucl -out $name" ; done

makeblastdb -in ../../../Data/Trinity-Build-Info/Individually/Mom/mom_trinity.Trinity.fasta -dbtype nucl -out mom

makeblastdb -in ../../../Data/Trinity-Build-Info/Individually/Dad/Dad_trinity.Trinity.fasta -dbtype nucl -out dad
```

2. Testing
```bash
blastn -task blastn -evalue 0.00005 -query Xl-genome-sox3.fa -db blast_dbs/mom -outfmt 6 > testing_blast_outs/mom_sox3.blasts

blastn -task blastn -evalue 0.00005 -query Xl-genome-sox3.fa -db blast_dbs/dad -outfmt 6 > testing_blast_outs/dad_sox3.blasts

for i in blast_dbs/BJE* ; do name=$(grep -o "BJE[0-9]*_[a-z]*" <(echo $i)) ; blastn -task blastn -evalue 0.00005 -query Xl-genome-sox3.fa -db blast_dbs/$name -outfmt 6 > testing_blast_outs/$name\_sox3.blasts ; done

# using my shitty blast parser for now to collect sequence
cd testing_blast_outs/
for i in BJE*blasts ; do name=$(grep -o "BJE[0-9]*_[a-z]*" <(echo $i)) ; blastParser.pl --blast=$i --fasta=../../../Data/Trinity-Build-Info/Individually/$name/$name\_trinity.Trinity.fasta --hit ; mv BlastSequences.fas $name\_sox3_blastHits.fa ; done

for i in BJE*fa ; do name=$(grep -o "BJE[0-9]*_[a-z]*" <(echo $i)) ; sed  -i "s/>/>$name/g" $i ; done
```

Hm, of course, blast is pulling out other genes (sox9) with similar domains. But, can't just take top trascriptome hit as there will be lots of top hits. Then again, could take the top hit of the gene set from each transcript and group together by similar top hits (should only be one from the genome). You know, a proper reciprocal blast approach.

**Blasting transcriptomes to the Xl chr8L sex linked gene set**

Need to set up a reciprocal blast, which will involve blasting each gene against the laevis genome. I guess I should use the whole genome just to trim down on noise (as outlined above). But, just grabbing the Xl genes from the genome, I lose positional information, so I'd have to fish that back out (with a best blast between the Xl genes of the sex linked region or whatever). Or, I need to make gene headers with positional information to grep it out. 

Each transcriptome to genome:
```bash
blastn -task blastn -evalue 0.0000000001 -db
```
