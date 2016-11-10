#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;
use Data::Dumper qw(Dumper);
use JSON;

my $seqname;
my $seq;
my $kmer;
my %kmerfreq;

my $ntIUPAC = "ATGCURYSWKMBDHVNatgcuryswkmbdhvn"; 
my $seqcounter = 0;

my $ksize = shift @ARGV;
my $file = shift @ARGV;
my $jsondb = shift @ARGV;

unless ($ksize && $file && $jsondb ) { die "You must specify a kmer size, and fasta file.\n"; }
my $filestream = Bio::SeqIO->new(
	-file => $file,
	-format => 'fasta'
	); 
	
open(JSONFILE,">",$jsondb) || die "Can't open file handle for $jsondb\n";
    
while( my $seqobj = $filestream->next_seq){
	my $time = localtime;
	$seqcounter++;
	$seq = $seqobj -> seq();
	$seqname = $seqobj->display_id();
	
	#regex approach - slower, higher overhead -  but more flexible
	#while ($seq =~ /(?=([ATGCNatgcn]{$ksize}))/g){ $kmercounts{ $1 }++; }

	#substring approach, faster leaner brute force approach
	my $i = 0;
	while (my $kmer =  substr($seq, $i, $ksize) ){
		if(length($kmer)<$ksize){ last; } #skip kmers at the burnt ends of sequences.
		if($kmer =~ m/[^$ntIUPAC]/){ die "Matched kmer $kmer with non-IUPAC nucleotide code in ",$seqname,"\n"; }
		my $revcomp = reverse_complement_IUPAC( $kmer);
		if(defined( $kmerfreq{ $revcomp })){ $kmerfreq{ $revcomp }{ $seqname }++; }
		else{ $kmerfreq{ $kmer }{ $seqname }++; }
		$i++; 
	}
 	print "$time : $seqcounter : ",scalar(keys(%kmerfreq))," ",$ksize,"mers were identified : $seqname\n";
	#if( $seqcounter > 5 ){ last; }
}

#print Dumper \%kmerfreq;

#print $seqcounter," sequences were indexed, ",scalar(keys(%kmercounts))," ",$ksize,"mers were identified\n";
#foreach my $key (sort {$kmercounts{$b} <=> $kmercounts{$a} } keys %kmercounts){ 
#	print $key,"\t",$kmercounts{ $key },"\n";
#}

#my $kmerjsondb = $json->encode( \%kmerfreq ); # print the HASH to JSON file
#print KMERDB encode_json( \%kmerfreq );
print JSONFILE encode_json( \%kmerfreq );
close JSONFILE;

sub reverse_complement_IUPAC {
        my $dna = shift;

        # reverse the DNA sequence
        my $revcomp = reverse($dna);

        # complement the reversed DNA sequence
        $revcomp =~ tr/ATGCURYSWKMBDHVNatgcuryswkmbdhvn/TACGAYRWSMKVHDBNtacgayrwsmkvhdbn/;
        return $revcomp;
}


