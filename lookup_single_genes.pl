#!/usr/bin/perl

use strict;
use warnings;
$ENV{BIOPERL_INDEX_TYPE} = "SDBM_File";
use Bio::Index::Fasta;
use Bio::SeqIO;

my $usage = '

Prior to runnin this script, make a master fasta db (ALL.fasta) by concatenating all dbs
used in the blast searches.

USAGE: perl lookup_single_genes.pl fastadb	

';

my $list = 'seq_list_toadd.txt'; #change this if it's not the correct name
my $fastadb = shift or die ("no database provided");
open IN, "<$list" or die "no $list file provided";

#make index of the fasta db, if already made skip it
my $idx = Bio::Index::Fasta->new(-filename => $fastadb . ".idx", -write_flag => 1);
$idx->make_index($fastadb);

while (my $line = <IN>) {
	chomp $line;
	my ($aln,$ac) = split (/\t/, $line);
	print "$aln\t$ac.....DONE\n";
	my $out = Bio::SeqIO->new (-format => 'Fasta', -file => ">>$aln");
	my $seq_obj = $idx->fetch($ac);
	$out->write_seq($seq_obj);
}

close IN;

exit;