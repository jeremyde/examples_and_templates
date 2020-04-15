#!/usr/bin/perl

#A script to bin markers in RIL or CSSLs.
#This is needed for dense genotyping
#Markers should be imputed first for best results.  Missing data creates bins.

#inputvariables tab delimited marker file and maximum distance of flanking markers to use for imputation

#marker file layout (tab delimited)
## marker individual1 individual2 etc.
#marker names as chr underscore position e.g., S1_10000
#markers must be coded as parent 1 (0) and partent 2 (2) and missing (-1)

#NOTE!! Bins will overlap because of recombination breakpoint uncertainty.  

use strict;
use Getopt::Std;
use Data::Dumper;

#reading options
our ($opt_i);
getopts('i:');
if (!$opt_i) {
    print STDERR "\nInput file (-i) required\n\n\n";
}

my @accessions;
my @markers;

my %markerchr;
my %markerpos;

my $firstline = 1;
my $secondline = 1;
my @prev_row;
my $bin_start_marker;
my $bin_end_marker;

open INFILE, "<", $opt_i or die "No such input file $opt_i";
while (<INFILE>) {
    chomp $_;
    if ($firstline == 1) {
	my @header = split('\t',$_);
	my @header = @header[1.. $#header];
	print "SNP\tRange\t".join("\t",@header)."\n";
	$firstline = 0;
	next;
    }
    my @row =  split('\t', $_);
    my $marker = $row[0];
    push @markers, $marker;
    my @markerinfo = split('_',$marker);
    $markerchr{$marker} = $markerinfo[0];
    $markerpos{$marker} = $markerinfo[1];
    @row = @row[1 .. $#row];

    if ($secondline == 1) {
	$bin_start_marker = $marker;
	$bin_end_marker = $marker;
	@prev_row = @row;
	$secondline = 0;
    }

    if (join(",",@row) eq join(",",@prev_row)) {
	$bin_end_marker = $marker;
	#note, this will misbehave if end of one chromosome has exact same marker scores as the start of the next
    }
    else {
	my @startmarkerinfo = split('_',$bin_start_marker);
	my @endmarkerinfo = split('_',$bin_end_marker);
	my $centerpos = ($markerinfo[1] - $startmarkerinfo[1])/2 + $startmarkerinfo[1];
	my $closest_to_center=10000000000;
	my $centermarker;
	foreach my $checkmarker (@markers) {
	    if ($markerchr{$checkmarker} eq $markerinfo[0]) {
		my $markerdist = abs($markerpos{$checkmarker} - $centerpos);
		    if ( $markerdist < $closest_to_center) {
			$closest_to_center = $markerdist;
			$centermarker = $checkmarker;
		}
	    }
	}
	print $centermarker."\t".$markerinfo[0]."_".$startmarkerinfo[1]."-".$markerinfo[1]."\t";
	print join("\t",@prev_row)."\n";
	@prev_row = @row;
	$bin_start_marker = $bin_end_marker;
	$bin_end_marker = $marker;
    }
}
close INFILE;

# deal with last line
my @startmarkerinfo = split('_',$bin_start_marker);
my @endmarkerinfo = split('_',$bin_end_marker);
my $centerpos = ($startmarkerinfo[1] - $endmarkerinfo[1])/2 + $startmarkerinfo[1];
my $closest_to_center=10000000000;
my $centermarker;
foreach my $checkmarker (@markers) {
    if ($markerchr{$checkmarker} eq $startmarkerinfo[0]) {
	my $markerdist = abs($markerpos{$checkmarker} - $centerpos);
	if ( $markerdist < $closest_to_center) {
	    $closest_to_center = $markerdist;
	    $centermarker = $checkmarker;
	}
    }
}
print $centermarker."\t".$startmarkerinfo[0]."_".$startmarkerinfo[1]."-".$endmarkerinfo[1]."\t";
print join("\t",@prev_row)."\n";

