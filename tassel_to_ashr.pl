#!/usr/bin/perl

use strict;
use Getopt::Std;
use Math::Complex;
use List::Util qw(sum);

sub mean {
    return sum(@_)/@_;
}


our ($opt_i, $opt_e);
getopts('i:e:');
if (!$opt_i) {
    print STDERR "\nInput file name (-i) required\n\n\n";
}
if (!$opt_e) {
    print STDERR "\nInput file name (-e) required\n\n\n";
}

open INFILE, "<", $opt_i or die "No such input file $opt_i";
open INFILE2, "<", $opt_e or die "No such input file $opt_e";

my $first_line = 1; #used to skip the header;
my $waiting_for_next_line = 0;
my $last_marker_seen;
my %effects_lookup;

# File 2 headers
# Trait   Marker  Locus   Site    Allele  Effect  Obs

print STDERR "Loading allele effects file...\n";
#load file 2 into hash
while (<INFILE2>) {

    my %marker_effect_info;
 
    chomp $_; #strips the new line character from the current input

    # turns the current tab delimited line into an array of values
    my @row =  split('\t', $_);

    #skips the header
    if ($first_line == 1) { 
	$first_line = 0;
	#check for the correct header
	if ($row[5] ne "Effect") {
	    die "Input file -e has the wrong header\n";
	}
	next;
    }

    my $marker_name = $row[1];

    if ($waiting_for_next_line == 0) {
	#starting a new marker entry
	
	$marker_name =  $row[1];
	if ($marker_name eq $last_marker_seen) {
	    die "Only bi-allelic markers allowed: $marker_name\n";
	}

	$marker_effect_info{'marker_allele_1'} =  $row[4];
	$marker_effect_info{'marker_effect'} =  $row[5];
	$marker_effect_info{'marker_obs_n1'} =  $row[6];

	if ($effects_lookup{$marker_name}) { 
	    die "Markers in wrong order in -e file at marker $marker_name\n";
	}
	#####
	$effects_lookup{$marker_name} = \%marker_effect_info;
	$waiting_for_next_line = 1;
	$last_marker_seen = $marker_name;
    } elsif ($last_marker_seen eq $marker_name) {  
	$waiting_for_next_line = 0;
	$last_marker_seen = 0;
	$effects_lookup{$marker_name}->{'marker_allele_0'} = $row[4];
	$effects_lookup{$marker_name}->{'marker_obs_n0'} = $row[6];
    } else {
	die "Did not find second allele for marker $marker_name\n";
    }
   
}
close INFILE2;
print STDERR "Loaded allele effects file\nLoading p-value file\n";


$first_line = 1;
my @errordf;

# Trait   Marker  Chr     Pos     df      F       p       add_effect      add_F   add_p   dom_effect      dom_F   dom_p   errordf MarkerR2        Genetic Var     Residual Var    -2LnLikelihood


print "Trait\tMarker\tChr\tPos\tp\tEffect\tStd_Err\t\Allele1\tAllele0\n";


while (<INFILE>) {

    chomp $_; #strips the new line character from the current input

    # turns the current tab delimited line into an array of values
    my @row =  split('\t', $_);

    #skips the header
    if ($first_line == 1) { 
	$first_line = 0;
	#check for the correct header
	if ($row[16] ne "Residual Var") {
	    die "Input file -i has the wrong header\n";
	}
	next;
    }

    my $marker_name = $row[1];

    my $pval = $row[6];
    if ($pval eq "NaN") { 
	next;
    }

    push @errordf, $row[13];

    my $marker_lookup = $effects_lookup{$marker_name};

    if ($marker_lookup) {
    
	print $row[0]."\t".$row[1]."\t".$row[2]."\t".$row[3]."\t".$row[6]."\t";
	my $marker_effect = $marker_lookup->{'marker_effect'};
	my $res_var = $row[16];
	my $stdev = sqrt($row[15] + $res_var);
	#my $zscore = $marker_effect/$stdev;
	#print $zscore."\t";
	
	my $obs_1 = $marker_lookup->{'marker_obs_n1'};
	my $obs_0 = $marker_lookup->{'marker_obs_n0'};
	my $std_err = sqrt($res_var/$obs_1 + $res_var/$obs_0);
	print $std_err."\t".$marker_lookup->{'marker_allele_1'}."\t".$marker_lookup->{'marker_allele_0'}."\n";
	

    } else {
	print STDERR "Warning:  Could not find allele effects for marker $marker_name\n";
    }
}	

print STDERR "Mean df of error: ".mean(@errordf)."\n";

close INFILE;

