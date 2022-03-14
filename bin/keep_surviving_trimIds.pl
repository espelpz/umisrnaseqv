#!/usr/bin/perl 

use strict;
use warnings;

my $x = $ARGV[0];
my $y = $ARGV[1];

# Get IDs from the reverse trim file

my %ids = ();

# open file with the trimmed reads to only keep the read IDS
#open(TRIMFILE, "zcat /mnt/lustre/groups/hodgkinsonlab/gcarbajo/projects/ctDNA/analysis/trimgalore/${sam}/*2_val_2.fq.gz|") or die "Cannot open UMI file $sam: $!\n";
open(TRIMFILE, "zcat $x |") or die "Cannot open trimmed reads file $x: $!\n";


my $counter1=0;
while(my $line1=<TRIMFILE>){
	$counter1++;
	if($counter1 == 1) {
		$ids{$line1}++;
		#print $line1;
		next;
	}	
	if ($counter1 == 4){
		$counter1 = 0;
	}
}

close TRIMFILE;

#open file of UMIs to only keep the one for which we have read IDS
#open(UMIFILE, "zcat /mnt/lustre/groups/hodgkinsonlab/gcarbajo/projects/ctDNA/data/ctDNA_IBIMA/HUMfeqzR/UMI/21060*/${sam}*fq.gz|") or die "Cannot open UMI file $sam: $!\n";
open(UMIFILE, "zcat $y |") or die "Cannot open UMI file $y: $!\n";

my $counter2=0;
while ( my $line=<UMIFILE> ) {
	if ($counter2 == 3) {
		print $line;
		$counter2 = 0;
		next;
	}
	if ($ids{$line}) {
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
