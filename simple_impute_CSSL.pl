#!/usr/bin/perl

#A simple script to impute markers in CSSLs based on flanking markers.
#This works for dense genotyping with sporadic missing data

#inputvariables tab delimited marker file and maximum distance of flanking markers to use for imputation

#marker file layout (tab delimited)
## marker individual1 individual2 etc.
#marker names as chr underscore position e.g., S1_10000
#markers must be coded as parent 1 (0) and partent 2 (2) and missing (-1)


use strict;
use Getopt::Std;
use Data::Dumper;

#reading options
our ($opt_i, $opt_d);
getopts('i:d:');
if (!$opt_i) {
    print STDERR "\nInput file (-i) required\n\n\n";
}
if (!$opt_d) {
    print STDERR "\nDistance for imputation based on flanking (-d) required\n\n\n";
}

my @accessions;
my @markers;

my $firstline = 1;
my $secondline = 1;
open INFILE, "<", $opt_i or die "No such input file $opt_i";
while (<INFILE>) {
    chomp $_;
    if ($firstline == 1) {
	print $_."\n";
	$firstline = 0;
	next;
    }
    my @row =  split('\t', $_);
    push @markers, $row[0];

    @row = @row[1 .. $#row];
    my $midx = 0;
    foreach my $score (@row) {
	if ($secondline == 1) {
	    my @markerarray;
	    push @accessions, \@markerarray;
	}
	push @{$accessions[$midx]}, $score;
	$midx++;
    }
    $secondline = 0;
}
close INFILE;

print STDERR scalar(@accessions)." accessions\n";
print STDERR scalar(@markers)." markers\n";
#print STDERR Dumper($accessions[$#accessions]);
#print STDERR Dumper($accessions[0]);

my %markerchr;
my %markerpos;

foreach my $marker (@markers) {
    my @markerarray = split('_',$marker);
    $markerchr{$marker} = $markerarray[0];
    $markerpos{$marker} = $markerarray[1];
}

my $midx = 0;
foreach my $marker (@markers) {
    print $marker;
    my $accidx = 0;
    foreach my $accref (@accessions) {
	my $accession = @{$accref};
	my $score = @{$accessions[$accidx]}[$midx];
	if ($score == -1) {
	    my $chrname = $markerchr{$marker};
	    my $pos = $markerpos{$marker};
	    
	    my $backward_marker_score = "notfound";
	    my $forward_marker_score = "notfound";
	    #searching forward and backward for flanking marker scores;
	    my $mcounter = 0;
	    while ($backward_marker_score eq "notfound") {
		$mcounter--;
		if ($midx+$mcounter == 0) { #don't search backwards on 1st marker
		    last;
		}
		
		my $flankmarker = $markers[$midx+$mcounter];
		my $flankchr = $markerchr{$flankmarker};
		my $flankpos = $markerpos{$flankmarker};
		my $flankscore = @{$accessions[$accidx]}[$midx+$mcounter];

		if ($chrname ne $flankchr) {
		    last;
		}
		if ($pos - $flankpos > $opt_d) {
		    last;
		}
		if ($flankscore ne "-1") {
		    $backward_marker_score = $flankscore;
		}
	    }
	    my $mcounter = 0;
	    while ($forward_marker_score eq "notfound") {
		$mcounter++;
		if ($midx+$mcounter == $#markers) { #don't search forwards on last marker
		    last;
		}
		
		my $flankmarker = $markers[$midx+$mcounter];
		my $flankchr = $markerchr{$flankmarker};
		my $flankpos = $markerpos{$flankmarker};
		my $flankscore = @{$accessions[$accidx]}[$midx+$mcounter];

		if ($chrname ne $flankchr) {
		    last;
		}
		if ($flankpos - $pos > $opt_d) {
		    last;
		}
		if ($flankscore ne "-1") {
		    $forward_marker_score = $flankscore;
		}
	    }

	    #test to see if flanking markers match score
	    unless ($forward_marker_score eq "notfound" || $backward_marker_score eq "notfound") {
		if ($forward_marker_score eq $backward_marker_score ) {
		    $score = $forward_marker_score;
		}
	    }
	}
	print "\t".$score;
	$accidx++;
    }
    print "\n";
    $midx++;
}

