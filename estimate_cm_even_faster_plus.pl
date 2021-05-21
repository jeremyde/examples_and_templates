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

my %mapmarkerchr;
my %mapmarkerpos;
my %mapmarkercm;


# push @markers,"S1_0";
# push @markers,"S2_0";
# push @markers,"S3_0";
# push @markers,"S4_0";
# push @markers,"S5_0";
# push @markers,"S6_0";
# push @markers,"S7_0";
# push @markers,"S8_0";
# push @markers,"S9_0";
# push @markers,"S10_0";
# push @markers,"S11_0";
# push @markers,"S12_0";


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
    #my @markerinfo = split('_',$marker);
    #my $chr = $markerinfo[0];
    #$chr =~ s/\D//g;
    #$markerchr{$marker} = $chr;
    #$markerpos{$marker} = $markerinfo[1];
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
    $mapmarkerchr{$marker} = $row[1];
    $mapmarkercm{$marker} = $row[2];
    $mapmarkerpos{$marker} = $row[3];

}
close MAPFILE;

print "Marker\tChr\tcM\tbp\n";
my $lastmapmarker = 0;
foreach my $marker (@markers) {
    my $startmarker;
    my $startmarker_pos;
    my $startmarker_cm;
    my $endmarker;
    my $endmarker_pos;
    my $endmarker_cm;
    my @markerinfo = split('_',$marker);
    my $chr = $markerinfo[0];
    $chr =~ s/\D//g;
    my $pos = $markerinfo[1];
    #$markerchr{$marker} = $chr;
    #$markerpos{$marker} = $markerinfo[1];

    for (my $i=$lastmapmarker; $i<=$#mapmarkers; $i++){
	my $mapmarker = $mapmarkers[$i];
	#print "c1:".$markerchr{$marker}."c2:".$markerchr{$mapmarker};
	if ($chr ne $mapmarkerchr{$mapmarker}) {
	    next;
	}
	my $mapmarkerpos_mapmarker = $mapmarkerpos{$mapmarker};
	if ($mapmarkerpos_mapmarker >= $pos) {
	    $endmarker = $mapmarker;
	    $endmarker_pos = $mapmarkerpos_mapmarker;
	    $endmarker_cm = $mapmarkercm{$endmarker};
	    $lastmapmarker = $i-1;
	    last;
	}
	$startmarker = $mapmarker;
	$startmarker_pos = $mapmarkerpos_mapmarker;
	$startmarker_cm = $mapmarkercm{$startmarker};
    }
    my $markercm;
    unless ($startmarker) {
	$markercm = $pos * ($endmarker_cm/$endmarker_pos);
    }
    unless ($endmarker) {
	$markercm = $pos * ($startmarker_cm/$startmarker_pos);
    }
    #print "st:$startmarker\n"; 
    unless ($markercm) {
	my $cm_per_bp = ($endmarker_cm - $startmarker_cm) / ($endmarker_pos - $startmarker_pos + 1);
	$markercm = $startmarker_cm + (($pos - $startmarker_pos) * $cm_per_bp);
    }
    if ($markercm < 0.0000000001) { $markercm = 0.0000000001;}  
    my $markercm_decimal = sprintf("%.10f", $markercm);
    print "$marker\t".$chr."\t".$markercm_decimal."\t".$pos."\n";
}
