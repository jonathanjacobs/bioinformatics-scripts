#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;

my $seqname;
my $seq;
my $kmer;
my @kmers;
my %kmercounts;

my $ksize = shift @ARGV;
my $file = shift @ARGV;
unless ($ksize && $file) { die "You must specify a kmer size and file name.\n"; }
my $filestream = Bio::SeqIO->new(
	-file => $file,
	-format => 'fasta'
	); 
 
my $seqcounter = 0;
while( my $seqobj = $filestream->next_seq){
	print "Counting ",$seqobj->display_id(),"\n";
	$seqcounter++;
	$seq = $seqobj -> seq();
	
	#regex approach - slower, higher overhead -  but more flexible
	#while ($seq =~ /(?=([ATGCNatgcn]{$ksize}))/g){ $kmercounts{ $1 }++; }
	
	#substring approach, faster leaner brute force approach
	my $i = 0;
	while (my $kmer =  substr($seq, $i, $ksize) ){
		if(length($kmer)<$ksize){ last; }
		$kmercounts{ $kmer }++;
		$i++; 
	}
	
	#array approach
	# TBD
	
	#last;
}

print $seqcounter," sequences were indexed, ",scalar(keys(%kmercounts))," ",$ksize,"mers were identified\n";
#foreach my $key (sort {$kmercounts{$b} <=> $kmercounts{$a} } keys %kmercounts){ 
#	print $key,"\t",$kmercounts{ $key },"\n";
#}
	
