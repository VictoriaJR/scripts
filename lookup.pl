#!/usr/bin/perl -w

##Script to look up a set of sequences from a fasta file
##reading from a list of keywords (id) in flat text.

##USAGE: perl lookup.pl fastafile id_list

use strict;
use warnings;
use Bio::SeqIO;

##Specify files on the command line. The $fasta is the file
##containing all sequences and $ids is the list of keywords/descriptions.  
my $fasta = shift;
my $ids = shift;

my $in  = Bio::SeqIO->new(-file => "$fasta", -format => 'Fasta');
my $out = Bio::SeqIO->new(-file => '>lookup_out.fasta', -format => 'Fasta');

my %hashids;

open IDS, "<$ids" or die "Can't open $ids : $!";

##Loop to populate hash from file IDS
while (my $line = <IDS>) {
	chomp $line;
	$hashids{$line} = 1;
}

close IDS;

##Go through each sequence in the fasta file
while (my $seq = $in->next_seq()) {
	my $id = $seq->id();
	if ($hashids{$id}) {
    	$out->write_seq($seq); 
	}
}

exit;
