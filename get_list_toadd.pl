#!/usr/bin/perl

use strict;
use warnings;

my $usage = '

This script is involved in the updating of the single-genes for phylogenomics.

It goes through the parsed blast output files, and pull out a list of sequence
accessions to add to the single-genes.

USAGE: perl get_list_toadd.pl *_parsed.txt

';

die $usage unless @ARGV;

open OUT, ">>seq_list_toadd.txt";

while (my $blast_parsed = shift @ARGV) {
	open IN, "<$blast_parsed" or die "can't open $blast_parsed";
	my @blast_file_name = split ('__', $blast_parsed); #split the long parsed file name to get only the name of the alignment file
	my $gene_name = $blast_file_name[0];
	print "Gene name: $gene_name\n";
	my %seen = ();
	while (my $line = <IN>) {
		chomp $line;
		my ($seq_id,$seq_length,$hit_id,$score,$evalue) = split (/\t/, $line);
		unless ($seen{$hit_id}) { #print OUT the different accession names only once
			if ($hit_id eq 'No hits found') { #do nothing if "no hit"
				next;
			}
			else {
				print OUT "$gene_name\t$hit_id\n";
				$seen{$hit_id}++;
			}
		}
	}
	close IN;
}

close OUT;

exit;