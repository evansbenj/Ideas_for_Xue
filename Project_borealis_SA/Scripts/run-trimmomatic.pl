#!/usr/bin/perl

use strict;
use warnings;

########
# Runs trimmomatic on paired files.
########

# my @IDs = ("BJE3896","BJE3897","BJE3929","BJE4009","BJE4017","BJE4039","BJE4072","BJE4082")

my $r1_name;
my $r2_name;
my $r1_par_out;
my $r1_unpar_out;
my $r2_par_out;
my $r2_unpar_out;
my $bash_IDs;
my @IDs;


my ($dir) = @ARGV;

$bash_IDs = `ls $dir | grep -o "BJE[0-9]*_[a-z]*_[a-z]*" | sort | uniq`;

@IDs = split("\n",$bash_IDs);

for my $name (@IDs) {
	$r1_name = $name."*R1.fastq.gz";
	$r2_name = $name."*R2.fastq.gz";
	$r1_par_out = $name."_R1_paired.fastq.gz";
	$r1_unpar_out = $name."_R1_unpaired.fastq.gz";
	$r2_par_out = $name."_R2_paired.fastq.gz";
	$r2_unpar_out = $name."_R2_unpaired.fastq.gz";

	system("java -jar /home/benf/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 $r1_name $r2_name $r1_par_out $r1_unpar_out $r2_par_out $r2_unpar_out ILLUMINACLIP:/home/benf/bin/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36");
}

exit;
