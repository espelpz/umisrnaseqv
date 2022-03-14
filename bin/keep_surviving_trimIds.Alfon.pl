#!/usr/bin/perl 

use strict;
use warnings;

my $x = $ARGV[0];
my $y = $ARGV[1];

# Get IDs from the reverse trim file

my %ids = ();

# open file with the trimmed reads to only keep the read IDS
#open(TRIMFILE, "zcat /mnt/lustre/groups/hodgkinsonlab/gcarbajo/projects/ctDNA/analysis/trimgalore/Alfon/${x}/*_val_2.fq.gz|") or die "Cannot open UMI file $x: $!\n";
open(TRIMFILE, "zcat $x |") or die "Cannot open UMI file $x: $!\n";
#open(TRIMFILE, "zcat  test_trimmed_reads.fq.gz|") or die "Cannot open UMI file $sam: $!\n";


my $counter1=0;
while(my $line1=<TRIMFILE>){
	$counter1++;

    my ($readId1 , $readId2) = split " " , $line1;

	if($counter1 == 1) {
    
		$ids{$readId1}++;
#		print $line1;
#    print "($readId1 , $readId2)\n";
#        print $readId1,"\n";
		next;
	}	
	if ($counter1 == 4){
		$counter1 = 0;
	}
}

close TRIMFILE;

#exit;

#open file of UMIs to only keep the ones for which we have read IDS
#open(UMIFILE, "zcat /mnt/lustre/groups/hodgkinsonlab/gcarbajo/projects/ctDNA/data/ctDNA_IBIMA/AlfonFastqfiles/${sam}*R2_001.fastq.gz|") or die "Cannot open UMI file $sam: $!\n";
open(UMIFILE, "zcat $y |") or die "Cannot open UMI file $y: $!\n";


my $counter2=0;
while ( my $line=<UMIFILE> ) {

    my ($readId1 , $readId2) = split " " , $line;    
#    print "($readId1 , $readId2)\n";
	if ($counter2 == 3) {
		print $line;
		$counter2 = 0;
		next;
	}
#	if ($ids{$line}) {
	if ($ids{$readId1}) {
		print $line;		
		$counter2++;
		next;
	}
	if ($counter2) {
		print $line;
		$counter2++;
		next;
	}
	
}	

close UMIFILE;

exit;
