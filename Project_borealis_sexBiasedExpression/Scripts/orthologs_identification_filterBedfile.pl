#!/usr/bin/perl

use strict;
use warnings;

my ($inputFile)=@ARGV;
my $blastout_file = "/home/xue/borealis_DE/de_sex_liver/mapping_blastn/liver_mvsf_blastout";
my %name; # structure: $name {transcript_id}
my %problem;
my @line;
my $num;


open(INPUT,"<", "$inputFile") or die "could not open the input file";
while (my $line = <INPUT>){
        chomp($line);
        @line = split (/\s+/,$line);
        
        #select the one with the highest mapping quality score; if with the same score, keep all of them 
        if (exists $name{$line[3]}){
        	if ($line[4] > $name{$line[3]}{"1"}[4]){ 
        		$name{$line[3]}{"1"}=[@line];
        	}
        	else{
        		$num = 1 + (scalar keys %{$name{$line[0]}});
        		$name{$line[3]}{"$num"}=[@line];
				if ($line[0] ne $name{$line[3]}{"1"}[0]){
        			$problem{$line[3]} = "different chr";
        		}
        		elsif($line[0] eq $name{$line[3]}{"1"}[0]){
        			$problem{$line[3]} = "same chr" unless (exists $problem{$line[3]} &&  $problem{$line[3]} ne "different chr");
        		}
        	}
        }
        else{
        	$name{$line[3]}{"1"}=[@line];
        }
}
close INPUT;

#mroe filtering:
#1. filter base on blast output
#2. filter base on genomic coordinate
my $num_alignment = 0;
my $counter_different =0;
my $counter_same=0;
my $best_hit_chr;
foreach my $id (sort {$a cmp $b} keys %name) {
	$num_alignment = (scalar keys %{$name{$id}});
	if ($num_alignment>1){
		$counter_different++;
	}

	if (exists $problem{$id} && $problem{$id} eq "different chr"){
		#$counter_different++;
		#if a transcript mapped to the different chr# with the same mapping score, check the blastn it against the blastn list
		open(INPUT,"<", "$blastout_file") or die "could not open the input file"; #check if a blastn output is available;
		READ: while (my $line = <INPUT>){
			if ($line[0] eq $id){
				$best_hit_chr = $line[1];
				last READ;
			}	
		}	
		close INPUT;

		foreach my $number (sort {$a cmp $b} keys %{$name{$id}}) {
			if ($name{$id}{$number}[0] ne $best_hit_chr){
				delete $name{$id}{$number};
			}
		}

	}

	#the alignments have the same mapping score and mapped to the same chr, keep the alignment with the lowest genomic coordinate (closer to the starting point)
	$num_alignment = scalar keys %{$name{$id}}; #check again after filtering base on blastn
	my $min_start =0;
	if (exists $problem{$id} && $problem{$id} eq "same chr" && $num_alignment>1){
		$counter_same ++;
		foreach my $number (sort {$a cmp $b} keys %{$name{$id}}) {
			if ($min_start && $min_start > $name{$id}{$number}[1]){
				$name{$id}{"1"} = [@{$name{$id}{$number}}];
				$min_start = $name{$id}{$number}[1];
			}
			else{
				$min_start = $name{$id}{$number}[1];
			}
		}
	}
}






#print the filtered alignment
PRINT: foreach my $id (sort {$a cmp $b} keys %name) {
	#print join ("\t", @{$name{$id}{"1"}}), "\n";
	foreach my $number (sort {$a cmp $b} keys %{$name{$id}}) {
	 	print join ("\t", @{$name{$id}{"1"}}), "\n";
	 	next PRINT;
	}
}

print "mapped to multiple chr: $counter_different; duplicates: $counter_same.\n";
