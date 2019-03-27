#!/usr/bin/perl

use strict;
use warnings;

my ($inputFile)=@ARGV;
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







#print the filtered alignment
PRINT: foreach my $id (sort {$a cmp $b} keys %name) {
	foreach my $number (sort {$a cmp $b} keys %{$name{$id}}) {
	 	print join ("\t", @{$name{$id}{"1"}}), "\n";
	 	next PRINT;
	}
}

print "mapped to multiple chr: $counter_different; duplicates: $counter_same.\n";
