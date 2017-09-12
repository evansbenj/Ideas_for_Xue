#!/usr/bin/perl

use strict;
use warnings;


##########
# Run bam adding read groups in a loop
##########

my ($dir, $who) = @ARGV; # who refers to some subset of the files in a dir, so that this runs on just them. will need to be unique, as it matches on a regex.

my $name;
my @Name;
my $RG;
my $left;
my $right;

opendir DH, $dir or die "could not open the directory\n";
foreach my $file (readdir DH) { # shit, has to handle R1 and R2
	if (($file =~ m/$who/) && ($file =~ /R1/)) {
		@Name = split("_",$file);
		$name = join("_",@Name[0..2]); # parse apart the file name

		$RG = "\'\@RG\\tLB:LB-$name\\tPL:ILLUMINA\\tID:CB2H9ANXX\\tSM:$name\'";

		$left = $file;
		$right = $name."_R2_scythe.fastq.gz";
		print "aligning $file\n";

		system("bwa mem -M -t 4 -R $RG /home/benf/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/Xla.v91.superscaffold.fa $dir/$left $dir/$right |  ~/bin/samtools-1.3.1/samtools view -b -F 4 - | ~/bin/samtools-1.3.1/samtools sort -o sort_$name.bam");
		# will do alignment, remove unmapped reads, sort the file, outputs bam.
	}
}
