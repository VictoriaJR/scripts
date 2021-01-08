#!/usr/bin/perl

use Bio::SearchIO;
use Bio::SeqIO;

use strict;
use warnings;

#Given a blast output in default format and the original fasta file, this script returns
#only the portion of the query that has a match.

#USAGE: perl rm_crap.pl fasta_file

my $fasta = $ARGV[0];
my $blastout = "$fasta".".sprot.blastout";

my $seqio = Bio::SeqIO->new(-file => $fasta, -format => "fasta");
my $blast_report = new Bio::SearchIO('-format' => 'blast', '-file' => $blastout);
open OUT, ">$fasta".".nocrap";	

my %sequences;
while (my $seq = $seqio->next_seq) {
	my $id = $seq->id;
	my $seq = $seq->seq();
	$sequences{$id} = $seq;
}	                     

#WARNING: set how many hits and hsps you want !

my $nbhits=1;
my $nbhsps=1;

while (my $result = $blast_report->next_result) {
	#output "no hits found" if that's the case
	if ($result -> num_hits == 0) {
		my $query_name = $result->query_name();
		print OUT ">$query_name\n$sequences{$query_name}\n";
   	}
   	else {  
		my $counter=0;
    	while (my $hit = $result->next_hit()) {
			if ($counter < $nbhits) {
				while (my $hsp = $hit->next_hsp()) {
					if ($counter < $nbhsps) {
						my $query_name = $result->query_name();
						my $start_query = $hsp->start('query');
						my $end_query = $hsp->end('query');
						my $length = (($end_query-$start_query)+1);
						my $short_seq = substr ($sequences{$query_name}, ($start_query-1), $length);
						print OUT ">$query_name\n$short_seq\n";
						print "$query_name\n";
						$counter+=1;
		  			}
				}
    		}
		}
	}	
}
exit;