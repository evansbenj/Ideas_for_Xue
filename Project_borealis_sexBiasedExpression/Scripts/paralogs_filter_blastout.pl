#!/usr/bin/perl

#################################################################Usage Information######################################################
#purpose: this script will filter the blast output to identify paralogs 
#usage: perl ~/script/paralogs
#approximate time cost: 

#input files: standard output from blastn output format6





#step used
#1.read the standard output, then store it 
#2.the script will go through the blast output, the first chromosome it matches to should be its own
###- if not match to its own, discard it. report no paralog found since the blastout does not find its match successfully
#3.the second chromosome or scaffold it matches to should be its paralogs. 
##- If the first match is correct and is on a chromosome, and it matches to its homologs, keep the paralog. 
##- If it is a scaffold, keep the second matches regardless of scasffold or chromosome. 
#4. go through the second match to identify the genomic location of the paralog
#5. output the location of its paralog
#
#another script, go through all the paralogs and identify the transcript bin that is within the genomic location of its paralog. 


########################################################################################################################################

use strict;
use warnings;

#my ($inputFile)=@ARGV;
my $blastout_file = "/home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_laevisGenomeApproach/laevis_gffgene_laevisV92_genome_blastout.tsv";
my %paralog; # structure: $gene{transcript_id}
my %temp;
my $chrName;
my $start;
my $end;
my $gene_id;
my  @line;
my $counter =0;
my %pGenelist;
my $max_diff=50; #the maximum of difference in query position to consider it as a different query
my $query_end_diff;
my $query_start_diff;
my $subject_start;
my $subject_end;


#filter the blast output (gff gene blastn to laevis genome) and keep the alignment that match to the homeologous chr or a scaffold  
open(INPUT,"<", "$blastout_file") or die "could not open the input file";
READ: while (my $line = <INPUT>){
      chomp($line);
      @line = split (/\s+/, $line);#split the line with spaces
      
      #if a paralog has been identified for the gene, compare 
      # if (exists $temp{$line[0]}){ 
      #           if ($line[1] eq $temp{$line[0]}[0]){
                      
      #                #print "$line\n";

      #                  $query_start_diff = abs ($temp{$line[0]}[3] - $line[6]); #the difference between the start and ending position of the query, if it is smaller than 10, it is likely the same section of the query
      #                  $query_end_diff = abs ($line[7] - $temp{$line[0]}[4]); #the difference between the start and ending position of the query, if it is smaller than 10, it is likely the same section of the query



      #                   if ($line[10] <= $temp{$line[0]}[5]){
      #                           if ($query_start_diff > $max_diff && $query_end_diff > $max_diff){
                                
      #                                   ($subject_start, $subject_end) = sort {$a <=> $b }($line[8], $line[9]); 

      #                                   if ($start < $temp{$line[0]}[1]){
      #                                           $temp{$line[0]}[1] = $subject_start;
      #                                   }

      #                                   if ($end > $temp{$line[0]}[2]){
      #                                           $temp{$line[0]}[2] = $subject_end;
      #                                   }   
      #                           }
      #                   }
      #           }
      #           else{
      #                   if ($paralog{$line[0]}){
      #                           next READ;
      #                   }
      #                   else{
      #                           $paralog{$line[0]}[0] = $temp{$line[0]}[0]; 
      #                           $paralog{$line[0]}[1] = $temp{$line[0]}[1];
      #                           $paralog{$line[0]}[2] = $temp{$line[0]}[2];
      #                           print "$line[0]\t$paralog{$line[0]}[0]\t$paralog{$line[0]}[1]\t$paralog{$line[0]}[2]\n";
      #                   }
      #           }
      #   #} 
      # }

      if ($paralog{$line[0]}){
                                next READ;
        }

      if (($gene_id) && ($line[0] eq $gene_id)){
        if ($line[1] eq $chrName){
                next READ;
        }
        else{
                #if the gene is from a scaffold, the second match should be its true paralog -> store it;
                #if the gene is from a chr, try to go through all the alignment until we find alignment that maps to its hommeologous chr;   
                if ($chrName =~ "caf"){ #the gene is from a scaffold, store the second seq match as paralog regardless it is mapped to a chr or another scaffold 
                        #$paralog{$gene_id}{"hemeoMatch"} = $line;  
                        ($start, $end) = sort {$a <=> $b }($line[8], $line[9]); 
                        $temp{$gene_id}[0] = $line[1];
                        $temp{$gene_id}[1]= $line[8];#start; 
                        $temp{$gene_id}[2] = $line[9];
                        $temp{$line[0]}[3] =$line[6];
                        $temp{$line[0]}[4] =$line[7];
                        $temp{$line[0]}[5] = $line[10]; #evalue

                        #############test#########
                                $paralog{$line[0]}[0] = $temp{$line[0]}[0]; 
                                $paralog{$line[0]}[1] = $temp{$line[0]}[1];
                                $paralog{$line[0]}[2] = $temp{$line[0]}[2];
                                print "$line[0]\t$paralog{$line[0]}[0]\t$paralog{$line[0]}[1]\t$paralog{$line[0]}[2]\n";
                        #############
                }
                elsif($chrName =~"hr" && $line[1] =~ "caf"){  #the gene is from a chromosome, alignment mapped to a scaffold
                        #$paralog{$gene_id}{"scaMatch"} = $line unless (exists $paralog{$gene_id}{"scaMatch"});
                }
                else{ #the gene is from a chromosome, alignment mapped to a chr
                        #print $gene_id, "\t", $line[1], "\t";
                        $line[1] = substr($line[1], 0, -1);
                        #print $line[1], "\n";
                        if ($chrName =~ $line[1]){
                                ($start, $end) = sort {$a <=> $b }($line[8], $line[9]); 
                                $temp{$gene_id}[0] = $line[1];
                                $temp{$gene_id}[1]= $start;#start; 
                                $temp{$gene_id}[2] = $end;
                                $temp{$line[0]}[3] =$line[6];
                                $temp{$line[0]}[4] =$line[7];
                                 $temp{$line[0]}[5] = $line[10]; #evalue
                                #print "$paralog{$gene_id}[0], $paralog{$gene_id}[1], $paralog{$gene_id}[2]\n";
                                #print $gene_id, "\t", $line[1], "\n";

                        #############test#########
                                $paralog{$line[0]}[0] = $temp{$line[0]}[0]; 
                                $paralog{$line[0]}[1] = $temp{$line[0]}[1];
                                $paralog{$line[0]}[2] = $temp{$line[0]}[2];
                                print "$line[0]\t$paralog{$line[0]}[0]\t$paralog{$line[0]}[1]\t$paralog{$line[0]}[2]\n";
                        #############
                        }
                        else{
                                #$paralog{$gene_id}{"otherChrMatch"} = $line unless(exists $paralog{$gene_id}{"otherChrMatch"});  
                        }
                }
        }
      }
      else{ #first alignment for the gene
        $line[0]=~ s/>//g;
        $gene_id = $line[0];
        ($chrName, $start,$end)=($gene_id=~/(.*?):(.*?)-(.*)/);# special edition for gff file seq

        ###########################put a check here to see if the first match correctly identify where the chr come from#################################
        if ($line[1] ne $chrName){
                #print "Problematic: $gene_id mapped to $line[1]";
                $counter ++;
        } 
      }
}

close INPUT;

exit;

#find out which gene the paralog alignment come from 
#read laevis gene gff file
my $inputFile = "/home/xue/genome_data/laevis_genome/gff3_xl9_2/XENLA_gene.gff3";
my $name;
my %gene;


open(INPUT,"<", "$inputFile") or die "could not open the input file";
while (my $line = <INPUT>){
        chomp($line);
        @line = split (/\s+/,$line);
        
        $name = "$line[0]:$line[3]-$line[4]";
       $gene{$line[0]}{$name}[0] = $line[3];
       $gene{$line[0]}{$name}[1] = $line[4];
}
close INPUT;
#print scalar keys %gene, "\n";


my %listbychr;
my %paraloglist;

foreach my $id (sort {$a cmp $b} keys %paralog) {
        # if (exists $paralog{$id}{"homeoMatch"}){
        #         #do nothing
        #         print $paralog{$id}{"homeoMatch"},"\n";
        # }
        # else{
        #         if (exists $paralog{$id}{"scaMatch"}){
        #                 #$paralog{$id}{"homeoMatch"} = $paralog{$id}{"scaMatch"};
        #                 #print $paralog{$id}{"homeoMatch"}, "\n";
        #         }
        #         else{
        #                 $paralog{$id}{"homeoMatch"} = "NA";
        #         }
        # }

         #print "$paralog{$id}\n";

       #@line = split (/\s+/, $paralog{$id}); 
      $chrName = $paralog{$id}[0];
      foreach my $name (sort {$a cmp $b} keys %{$gene{$chrName}}) {
                #if ($gene{$line[0]}{$name}[0] < $line[8] && $gene{$line[0]}{$name}[1] > $line[9]){
                if (($gene{$chrName}{$name}[0] < $paralog{$id}[1] && $gene{$chrName}{$name}[1] > $paralog{$id}[1]) ||
                ($gene{$chrName}{$name}[0] < $paralog{$id}[2] && $gene{$chrName}{$name}[1] > $paralog{$id}[2])){
                        $paraloglist{$id} = $name;
                        push (@{$listbychr{$line[1]}}, $name);
                        #print "$id\t$name\n";
                }
        }
}


#identify borealis transcript that is in the paralog region
$inputFile = "/home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_laevisGenomeApproach/borealis_denovoT_laevisV92_genome_gmap.bed";
open(INPUT,"<", "$inputFile") or die "could not open the input file";
while (my $line = <INPUT>){
        chomp($line);
        @line = split (/\s+/,$line);
        
        foreach my $name (@{$listbychr{$line[1]}}) {
                
                my ($chrname, $gene_start,$gene_end)=($name=~/(.*?):(.*?)-(.*)/);
                
                if ($gene_start < $line[1] && $gene_end > $line[2]){
                        push (@{$pGenelist{$name}{"borealis"}}, "xb_$line[3]");
                }
        }
}
close INPUT;


#identify laevis transcript that is in the paralog region
$inputFile = "/home/xue/borealis_DE/de_sex_liver/borealis_laevis_tropicalis_orthologs/orthologs_laevisGenomeApproach/laevis_denovoT_laevisV92_genome_gmap_bedfile_filtered.bed";
open(INPUT,"<", "$inputFile") or die "could not open the input file";
while (my $line = <INPUT>){
        chomp($line);
        @line = split (/\s+/,$line);
        
        foreach my $name (@{$listbychr{$line[1]}}) {
                
                my ($chrname, $gene_start,$gene_end)=($name=~/(.*?):(.*?)-(.*)/);
                
                if ($gene_start < $line[1] && $gene_end > $line[2]){
                        push (@{$pGenelist{$name}{"laevis"}}, "xl_$line[3]");
                }
        }
}
close INPUT;

#print the list of orthologs transcript ID
foreach my $gene (sort {$a cmp $b} keys %paraloglist) {
        my $paralog_gene = $paraloglist{$gene};
        print "$gene\t$paralog_gene\t";
        
        if (exists $pGenelist{$paralog_gene} && exists $pGenelist{$paralog_gene}{"borealis"}){
                print join (",", @{$pGenelist{$paralog_gene}{"borealis"}});
                print "\t"; 
        }
        else{
                print "NA\t";
        }

        if (exists $pGenelist{$paralog_gene} && exists $pGenelist{$paralog_gene}{"laevis"}){
                print join (",", @{$pGenelist{$paralog_gene}{"laevis"}});
                print "\n"; 
        }
        else {
                print "NA\n";
        }

}

exit;
