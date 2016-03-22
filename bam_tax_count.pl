use strict;
use warnings;
use Bio::DB::Taxonomy;
use Bio::DB::Taxonomy::flatfile;

# bam_tax_count.pl
# v. 1.1
# See POD documentation at the end -- or type "perldoc bam_tax_count.pl" at the CL

# reconfigure depending on your system
my $email = 'jonathan.jacobs@gmail.com';
my $taxonomydir = '/home/share/NCBI/taxonomy/';
my $nodesfile = $taxonomydir.'nodes.dmp';
my $namesfile = $taxonomydir.'names.dmp';

my %reads_to_gi; #had of read names mapped to a REF of GI's (which is an array)
my %gi_counts; #hash of Genome ID's => #mapped reads
my %taxon_counts; # hash of TAXONOMIC ID's => #mapped reads

# cant seem to get local flatfile to work... so, using Entrez method (MUCH SLOWER)
my $taxdb = Bio::DB::Taxonomy->new( -email	=> $email, -source => 'entrez', -location => 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/' ); 

#create hash of arrays of $reads_to_gi{ $READNAME } => [ GI1, GI2, GI3... GIn]
foreach my $bamline (<STDIN>) {
	my @cols = split("\t",$bamline);
	my ($gi) = $cols[2] =~ m/gi\|(\S+?)\|\S+?\|\S+?\|/g;
	push( @{ $reads_to_gi { $cols[0] } }, $gi ); #push the GI# onto the array of GI's for that specific read
}

#now count all reads for each GI. Multi-mapped reads get fractional values - ie if a read maps to 2 GIs, it gets 0.5 per GI
foreach my $read (keys %reads_to_gi){
	foreach my $gi (@{$reads_to_gi{$read}}){
		if( exists $gi_counts{ $gi } ){
			$gi_counts{ $gi } = $gi_counts{ $gi } + (1 / scalar(@{$reads_to_gi{$read}}));
		}else{
			$gi_counts{ $gi } = (1 / scalar(@{$reads_to_gi{$read}}))
		} ;
	}
}

# call rollup_taxons and add up the read counts
foreach my $gi (keys %gi_counts){
	my $taxon = $taxdb->get_taxon( -gi=> $gi, -db=> 'nucleotide' );
    rollup_taxons( $taxon, $gi_counts{ $gi });
}

#export the results
foreach my $taxonid (keys %taxon_counts ){
	my $taxon = $taxdb->get_taxon(-taxonid => $taxonid);
	my @names = @{$taxon->name('scientific')};
	
	print $taxonid,"\t";
	print $taxon->rank(),"\t";
	print $names[0],"\t";
	print $taxon_counts{ $taxonid },"\n";
}

exit;

#recursively roll up taxon counts
sub rollup_taxons{
	my ($taxon,$reads) = @_;

	if (exists($taxon_counts{ $taxon->id() })){
		$taxon_counts{ $taxon->id() } = $taxon_counts{ $taxon->id() } + $reads;
	}else{
		$taxon_counts{ $taxon->id() } = $reads;
	}
	
	if( defined( $taxon->parent_id() ) ){
		#print $taxon->parent_id(),"\t",$taxon->id(),"\n";
		my $parent_taxon = $taxdb->get_taxon( -taxonid=>$taxon->parent_id());
		rollup_taxons($parent_taxon,$reads);
	}
	
	return 1;
	
}

=pod

=head1 NAME

bam_tax_count.pl

=head1 DESCRIPTION

This script currently assumes "samtools view [bamfile]" as (piped) STDIN input

For example:

> samtools view myfile.bam | perl bam_tax_count.pl

Output is a tab delimited file of 

TaxonID	TaxonRank	TaxonName	#ReadsMapped

This script prints out all reads for each taxonomic level, and rolls them up to each successively higher taxon. Multimapped reads are equally split between multiple taxons, i.e. if 2 reads are mapped to two species, each species will get 0.5 reads mapped. This is shorthand method for handling multimapped reads - future implementations might consider weighting the fractional mapping count to alignment quality, etc.

=head1 REQUIREMENTS

=over 12

=item 1. Bio::DB::Taxonomy

=item 2. Bio::DB::Taxonomy::flatfile

=back
 
Internet access to NCBI is requred since this script uses the current ENTREZ taxonomy. Local flatfile calling of taxonomy is not currently supported.

=head1 AUTHOR

Jonathan Jacobs / jonathan.jacobs@gmail.com


=head1 LICENSE

FreeBSD License

Copyright (c) 2016, Jonathan Jacobs
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

=head1 CHANGELOG

CHANGE LOG

v1.1 (22 MAR 2016)
- added POD documentation

v1.0 (18 MAR 2016)
- initial version
- cleaned up code a bit 

=cut

