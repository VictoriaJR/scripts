#!/usr/bin/perl

### This script performs blast searches of a fasta file against several databases provided by the
### user.

### blastall needs to be in the PATH, and the path to the formatdb databases has to be correct.

### USAGE: perl running_blast.pl fasta_file blastdb_name

### fasta_file contains the file that we want to blast.
### blastdb_name is a file with the formatdb names of the databases (as they appear in the blastdb 
### folder), one on each line.

use strict;
use warnings;

my $fasta_file = shift;
my $blastdb_name = shift;

open (IN, "<$blastdb_name") or die ("cannot open $ARGV[1]");

while (my $blastdb=<IN>) {
	chomp $blastdb;

my $cmd = "blastp -query $fasta_file -num_threads 16 -evalue 1e-10 -num_descriptions 4 -num_alignments 4 -db /Data/Lily/Colp12_Project/Phylogenomics_Aug2016/databases_forBLAST/$blastdb -out $fasta_file"."__"."$blastdb.blastout";
	my $oops = system $cmd; #system call and save return value in $oops
	die "FAILED $!" if $oops;
}

exit;
