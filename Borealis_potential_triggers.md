# Using RNAseq from female and male tads, and also the genome seqs from a male and female individual, I identified transcripts that are expressed only in the female gonad and are also only found in the female genome.  These transcripts are not expressed in adult female liver, or in male tad gonads, male liver, or testis and they are not in the male genome.

I extracted the sequences from the trinity assembly like this, where 'ids.file` has the list of transcripts.

```
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' ids.file fasta.file
```
