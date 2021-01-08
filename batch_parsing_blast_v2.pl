#!/usr/bin/perl

### This script performs parsing of blast searches, several at a time if more than one
### is provided.

### USAGE: perl batch_parsing_blast_v1.pl *.blastoutput
### Only the best hit will be reported

#####################################################################################################
### UPDATE: this version takes into account the evalue (1e-20) and query and hit coverage (both >0.5)
#####################################################################################################

use strict;
use warnings;

use Bio::SearchIO;

my $blast_output;

my $usage = '

USAGE: perl batch_parsing_blast_v1 *.blastoutput

';

die $usage unless @ARGV;

my $nbhits = 1;
#my $nbhsps = 1;

######## Modify this for different cutoff ########
my $evalue_cutoff = "1e-20";
my $coverage_cutoff = "0.5";
##################################################

###$evalue <= $evalue_threshold

while (my $blast_output = shift @ARGV) {
	chomp $blast_output;
	my $report_object = Bio::SearchIO -> new (-format => 'blast', 
                          				  	  -file   => $blast_output);
	open (OUT, ">$blast_output"."_parsed.txt");
	#parse the blast output
	while (my $result_object = $report_object -> next_result) {
		my $query_name = $result_object -> query_name;
		my $query_length = $result_object -> query_length;
		my $hit_counter = 0;
		while (my $hit_object = $result_object -> next_hit and ($hit_counter < $nbhits)) {
			$hit_counter = $hit_counter+1;
			my $hit_name = $hit_object -> name;
			my $score = $hit_object -> raw_score;
			my $evalue = $hit_object -> significance;
			my $query_coverage = $hit_object->frac_aligned_query();
			my $hit_coverage = $hit_object->frac_aligned_hit();
			print "$query_name\t$evalue\t$query_coverage\t$hit_coverage\n";
			if (($evalue <= $evalue_cutoff) && ($query_coverage >= $coverage_cutoff) && ($hit_coverage >= $coverage_cutoff)) {
				print OUT "$query_name\t$query_length\t$hit_name\t$score\t$evalue\n";
			}
		}
	}
	close (OUT);
}
exit;