#!/usr/bin/perl

#inputvariables segment file made by gwas_segment_maker.pl and a "3.txt" allele effects file output from Tassel

#segmentfile header
#CHR     Start   Stop    Peak    Peak_val        Merged

#Tassel "3.txt" allele effects file
#Trait   Marker  Locus   Site    Allele  Effect  Obs

use strict;
use Getopt::Std;
use Data::Dumper;

#reading options
our ($opt_e, $opt_s);
getopts('e:s:');
if (!$opt_e) {
    print STDERR "\nAllele effects file (-e) required\n\n\n";
    die;
}

if (!$opt_s) {
    print STDERR "\Segment file (-s) required\n\n\n";
    die;
}

#open the allele effects file or die
open EFFECTSINFILE, "<", $opt_e or die "No such input file $opt_e";

my $seen_SNPs;
my $effect_lookup; #this is a pointer to the hash to lookup allele effects;
my $first_line = 1; #used to recognize the header;
#reads the effects file line by line
while (<EFFECTSINFILE>) {
    # the variable $_ contains the current line.
    chomp $_; #strips the new line character from the current input
    #skip the header row
    if ($first_line==1) {
	#TODO add some checks to make sure the header is what we expect
	$first_line=0;
	next;
    }
    my @row =  split('\t', $_);
    my $chromosome = $row[2];
    my $position = $row[3];

    unless($effect_lookup->{$chromosome}->{$position}) {
	$effect_lookup->{$chromosome}->{$position} = [];
    }

    push @{$effect_lookup->{$chromosome}->{$position}}, $_;

}
#print Dumper($effect_lookup);
close EFFECTSINFILE;

#open the segment file or die
open SEGMENTINFILE, "<", $opt_s or die "No such input file $opt_e";

#my $seen_SNPs;

$first_line = 1; 
#reads the segment file line by line
while (<SEGMENTINFILE>) {
    chomp $_; 
    if ($first_line==1) {
	#TODO add some checks to make sure the header is what we expect
	print $_."\tAllele1\tAllele2\tAllele1Count\tAllele2Count\tAllele1Effect\tAllele2Effect\tTotalAlleles\n";
	$first_line=0;
	next;
    }
    my @row =  split('\t', $_);
    my $chromosome = $row[0];
    my $position = $row[3];

    my @effectlines = @{$effect_lookup->{$chromosome}->{$position}};

    #check for multi-allelic SNPs
    my $total_allele_count = scalar(@effectlines)/2;
    
    my %snp_info;

    my $which_of_pair = 0;
    foreach my $effectline (@effectlines) {
	my @effect_row =  split('\t', $effectline);
	my $allele = $effect_row[4];
	my $effect = $effect_row[5];
	my $allele_count = $effect_row[6];

	
	if ($which_of_pair == 0) {
	    $snp_info{'1_allele'} = $allele;
	    $snp_info{'1_allele_effect'} = $effect;
	    $snp_info{'1_count'} = $allele_count;
	    $which_of_pair++;
	    
	} elsif ($which_of_pair == 1) {
	    $snp_info{'2_allele'} = $allele;
	    $snp_info{'2_allele_effect'} = $effect;
	    $snp_info{'2_count'} = $allele_count;
	    $which_of_pair = 0;
	    print $_."\t".
		$snp_info{'1_allele'}."\t".
		$snp_info{'2_allele'}."\t".
		$snp_info{'1_count'}."\t".
		$snp_info{'2_count'}."\t".
		$snp_info{'1_allele_effect'}."\t".
		$snp_info{'2_allele_effect'}."\t".
		$total_allele_count."\n";

	    
	} else {
	    die "ERROR: SNPs pairs out of sync\n";
	}
    }
    
}
