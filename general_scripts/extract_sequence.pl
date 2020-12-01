#!/usr/bin/perl
#
use strict;
use warnings; 
use Scalar::Util qw(looks_like_number);

#######################################################
#purpose: extract sequence from fasta file
#usage: perl extract_sequence.pl seq_name.file fasta.file position_of_sequence_in_seq_name.file > output
#
#
#####################################################



my ($inputfile1, $inputfile2, $ID_position) = @ARGV;#$ID_position is the # of the column that contains the transcript ID   
my %transcript = ();
my @line = ();

open(INPUT,"<", "$inputfile1") or die "could not open the input file1";#opent the first input file
#open INPUT,"samtools view $inputfile1 |";

$ID_position = $ID_position - 1; #-1 because perl start from position 0
while (my $line = <INPUT>){
        chomp($line);
        #@line = split (/\t/, $line);#split the line with spaces
   	@line = split (/,/, $line);#split the line with comma
	transcript{$line[$ID_position]}= 1;
}
close INPUT;



open(INPUT,"<", "$inputfile2") or die "could not open the input file1";#opent the first input file
my $print = 0;
while (my $line = <INPUT>){
      chomp($line);
      @line = split (/\s+/, $line);#split the line with spaces
      if ($line =~ ">"){ # check if it is the seq name`
   	$print =0;
	$line[0]=~ s/>//g; #extract the transcript name, need it for checking in if
        
	#special edition for gff file seq
	#my ($chr, $start,$end)=($line[0]=~/(.*?):(.*?)-(.*)/);
	#$start = $start +1;
	#$line[0] = "$chr:$start-$end";
	
	if (exists $transcript{$line[0]}){
          print "$line\n"; #print sequence name
          $print = 1; #set the printing of sequence to 1 so that it goes to elsif to print the seq 
        }
      }
      elsif($print == 1){
        print "$line\n"; #print seq
      }
}

close INPUT;

