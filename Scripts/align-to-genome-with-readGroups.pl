#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;

##########
# Run bam adding read groups in a loop
##########

my ($dir, $who) = @ARGV; # who refers to some subset of the files in a dir, so that this runs on just them. will need to be unique, as it matches on a regex.

my $name;
my @Name;
my $RG;
my $left;
my $right;
my $use_bwa;
my $use_hisat;
my $rg_id;
my $sm;
my $pl;
my $lb;
my $use_star;

GetOptions(
	'bwa' => \$use_bwa,
	'hisat' => \$use_hisat,
	'star' => \$use_star
);

opendir DH, $dir or die "could not open the directory\n";
foreach my $file (readdir DH) { # shit, has to handle R1 and R2
	if (($file =~ m/$who/) && ($file =~ /R1/)) {
		@Name = split("_",$file);
		$name = join("_",@Name[0..2]); # parse apart the file name

		$RG = "\'\@RG\\tLB:LB-$name\\tPL:ILLUMINA\\tID:CB2H9ANXX\\tSM:$name\'";

		# alt handling of RG info
		$rg_id = "CB2H9ANXX";
		$sm = "SM:$name";
		$pl = "PL:ILLUMINA";
		$lb = "LB:LB-$name";

		$left = $file;
		$right = $name."_R2_scythe.fastq.gz";
		print "aligning $file\n";

		# will do alignment, remove unmapped reads, sort the file, outputs bam.
		if ($use_bwa) {
			system("bwa mem -M -t 4 -R $RG /home/benf/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/Xla.v91.superscaffold.fa $dir/$left $dir/$right |  ~/bin/samtools-1.3.1/samtools view -b -F 4 - | ~/bin/samtools-1.3.1/samtools sort -o sort_$name.bam");
		}
		elsif ($use_hisat) {
			system("~/bin/hisat2-2.1.0/hisat2 -p 4 --rg-id $rg_id --rg $sm --rg $pl --rg $lb -x /home/benf/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/HISAT-Index/hisat-ref -1 $dir/$left -2 $dir/$right | ~/bin/samtools-1.3.1/samtools view -b -F 4 - | ~/bin/samtools-1.3.1/samtools sort -o sort_$name.bam ")
		}
		elsif ($use_star) {
			system("~/bin/STAR-2.5.3a/bin/Linux_x86_64/STAR --runThreadN 5 --genomeDir /home/benf/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/STAR-Index/ --readFilesCommand zcat --outFileNamePrefix $name-star --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 20 --outFilterMismatchNoverLmax 0.5 --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 --outSAMattrRGline ID:$rg_id $sm $pl $lb --readFilesIn $dir/$left $dir/$right")
		}
	}
}


# --runThreadN 10 --genomeDir /home/benf/Borealis_Genome_HiSeqX/Analyses/Laevis_SuperScaffold_Reference/STAR-Index/ --readFilesCommand zcat --outFileNamePrefix dad-star-lessStringent --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 20 --outFilterMismatchNoverLmax 0.5 --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 --readFilesIn ../Trimmed/BJE3896_dad_liver_R1_scythe.fastq.gz ../Trimmed/BJE3896_dad_liver_R2_scythe.fastq.gz

# --outSAMattrRGline ID:XXXXX PG:XXXXX etc.
