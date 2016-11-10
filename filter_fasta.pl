#!/usr/bin/perl -w
# 
# This script filters FASTA files for sequences with specific designations, patterns, etc and outputs to standard out.
#
# for example
#
# perl filter_fasta.pl bigfasta.fna NC_ > subset.fna
#
# would write all the files with NC_ in the name line to subset.fna 

use strict;
use Bio::SeqIO;

my $file = shift @ARGV;
my $pattern = shift @ARGV;
unless ($file && $pattern) { die "You must specify a file name and an accession name pattern to match.\n"; }

my $filestream = Bio::SeqIO->new( -file => $file, -format => 'fasta' ); 
 
while( my $seqobj = $filestream->next_seq){
	if($seqobj->display_id() =~ m/$pattern/){
		print ">",$seqobj->display_id(),"\n";
		print $seqobj->seq(),"\n";
		}
}

# TO DO
# rewrite more portable version where Bio::Perl is not available
