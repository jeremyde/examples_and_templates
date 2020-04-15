#!/usr/bin/perl

#A script to estimate cM from a map based on physicial position of markers

#inputvariables tab delimited marker file and map file

#marker file only uses first column for marker names

#map file layout (tab delimited)
## Marker Chr cM bp

use strict;
use Getopt::Std;
use Data::Dumper;

#reading options
our ($opt_i, $opt_m);
getopts('i:m:');
if (!$opt_i) {
    print STDERR "\nInput marker file (-i) required\n\n\n";
}
if (!$opt_m) {
    print STDERR "\nInput map file (-m) required\n\n\n";
}

my @markers;
my @mapmarkers;

my %markerchr;
my %markerpos;
my %markercm;

my $firstline = 1;

open INFILE, "<", $opt_i or die "No such input file $opt_i";
while (<INFILE>) {
    chomp $_;
    if ($firstline == 1) {
	$firstline = 0;
	next;
    }
    my @row =  split('\t', $_);
    my $marker = $row[0];
    push @markers, $marker;
    my @markerinfo = split('_',$marker);
    my $chr = $markerinfo[0];
    $chr =~ s/\D//g;
    $markerchr{$marker} = $chr;
    $markerpos{$marker} = $markerinfo[1];
}
close INFILE;

$firstline = 1;
open MAPFILE, "<", $opt_m or die "No such input file $opt_i";
while (<MAPFILE>) {
    chomp $_;
    if ($firstline == 1) {
	$firstline = 0;
	next;
    }
    my @row =  split('\t', $_);
    my $marker = $row[0];
    push @mapmarkers, $marker;
    $markerchr{$marker} = $row[1];
    $markercm{$marker} = $row[2];
    $markerpos{$marker} = $row[3];
}
close MAPFILE;

print "Marker\tChr\tcM\tbp\n";
foreach my $marker (@markers) {
    my $startmarker;
    my $endmarker;
    foreach my $mapmarker (@mapmarkers) {
	#print "c1:".$markerchr{$marker}."c2:".$markerchr{$mapmarker};
	if ($markerchr{$marker} ne $markerchr{$mapmarker}) {
	    next;
	}
	if ($markerpos{$mapmarker} >= $markerpos{$marker}) {
	    $endmarker = $mapmarker;
	    last;
	}
	$startmarker = $mapmarker;
    }
    my $markercm;
    unless ($startmarker) {
	$markercm = $markerpos{$marker} * ($markercm{$endmarker}/$markerpos{$endmarker});
    }
    unless ($endmarker) {
	$markercm = $markerpos{$marker} * ($markercm{$startmarker}/$markerpos{$startmarker});
    }
    #print "st:$startmarker\n"; 
    unless ($markercm) {
	my $cm_per_bp = ($markercm{$endmarker} - $markercm{$startmarker}) / ($markerpos{$endmarker} - $markerpos{$startmarker} + 1);
	$markercm = $markercm{$startmarker} + (($markerpos{$marker} - $markerpos{$startmarker}) * $cm_per_bp);
    }
    print "$marker\t".$markerchr{$marker}."\t".$markercm."\t".$markerpos{$marker}."\n";
}
