#!/usr/bin/perl

#inputvariables tab delimited list of files and experiment names

use strict;
use Getopt::Std;
use DateTime qw( );

#reading options
our ($opt_i);
getopts('i:');
if (!$opt_i) {
    print STDERR "\nInput file (-i) required\n\n\n";
}

my @files;
my @exp_names;
open INFILE, "<", $opt_i or die "No such input file $opt_i";
while (<INFILE>) {
    chomp $_;
    my @row =  split('\t', $_);
    push @files, $row[0];
    push @exp_names, $row[1];
}
close INFILE;

my @segments;
my $name_index = 0;
foreach my $file (@files) {
    #print STDERR "segfile $file\n";
    open SEGFILE, "<", $file or die "No such input file $file";
    my $first_line = 1;
    my $last_stop;
    my $last_chr;
    while (<SEGFILE>) {
	#skips the header
	if ($first_line) { 
	    $first_line = 0;
	    next;
	}
	chomp $_;
	my @row =  split('\t', $_);
	my %segment;
	$segment{'chr'} = $row[0];
	#todo: strip non numeric from chromosome names to enable correct sorting;
	$segment{'start'} = $row[1];
	$segment{'stop'} = $row[2];
	$segment{'peak'} = $row[3];
	$segment{'name'} = $exp_names[$name_index];
	#checking for correct order and no overlaps withing a trait
	if ($last_stop) {
	    if ($last_chr ne $segment{'chr'}) {
		$last_stop = 0;
	    } else {
		if ($last_stop > $segment{'start'} ) {
		    print STDERR "Error in input file $file\n Ranges should be in order and not overlap within an experiment\n";
		    die;
		}
	    }
	}
	$last_stop = $segment{'stop'};
	$last_chr = $segment{'chr'};
	push @segments, \%segment;
    }
    $name_index++;
    close SEGFILE;
}

#sort all segments by chromosome then start position
my @segments_sorted  = sort { $a->{'chr'} <=> $b->{'chr'} or $a->{'start'} <=> $b->{'start'} } @segments;

my $chromosome;
my $overlap_start;
my $overlap_stop;
my $overlap_group_info;
my $overlap_segments;

sub commify {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text
}

my $overlap_counter = 0;
sub print_overlaps {
    if ($overlap_group_info) {
	#only print overlap group if it has more than one member
	if (scalar @{$overlap_group_info->{'segments'}} > 1) {
	    #print stderr "Printing group\n";
	    $overlap_counter++;
	    print "\n\tChr\tStart\tEnd\t\Peak_SNP\tLength\tTrait\n";
	    
	    my @print_segments_sorted  = sort { $a->{'name'} cmp $b->{'name'} } @{$overlap_group_info->{'segments'}};
	    #foreach my $r_segment (@{$overlap_group_info->{'segments'}}) {
	    foreach my $r_segment (@print_segments_sorted) {
		print "\t".$r_segment->{'chr'}."\t".
		    commify($r_segment->{'start'})."\t".
		    commify($r_segment->{'stop'})."\t".
		    commify($r_segment->{'peak'})."\t".
		    commify($r_segment->{'stop'}-$r_segment->{'start'})."\t".$r_segment->{'name'}."\n";
	    }
	    print "Range:\t".$overlap_group_info->{'chr'}."\t".
		commify($overlap_group_info->{'start'})."\t".
		commify($overlap_group_info->{'stop'})."\t"."\t".
		commify($overlap_group_info->{'stop'}-$overlap_group_info->{'start'})."\t\n";
	}
    }
    $overlap_group_info = 0;
}

print "#Overlapping segments\n";
print "#Report for segments from files listed in $opt_i\n";
my $dt = DateTime->now( time_zone => 'local' );
print "#Date: ".$dt->strftime("%m/%d/%Y")."\n";

foreach my $segment (@segments_sorted) {
    if ($chromosome) {
	if ($chromosome ne $segment->{'chr'}) {
	    if ($overlap_group_info) { 
		print_overlaps();		
		$overlap_group_info = 0;
		$overlap_start = 0;
		$overlap_stop = 0;
	    } 
	}
    }
    $chromosome = $segment->{'chr'};
    if ($overlap_group_info) {
	if ($segment->{'start'} <= $overlap_group_info->{'stop'}) {
	    if ($overlap_group_info->{'stop'} < $segment->{'stop'}) {
		$overlap_group_info->{'stop'} = $segment->{'stop'};
	    }
	    push @{$overlap_group_info->{'segments'}}, $segment;
	    next;
	} else {
	    print_overlaps();
	}
    }
    my %ginfo;
    my @overlap_array;
    $overlap_group_info = \%ginfo;
    $overlap_group_info->{'chr'} = $segment->{'chr'};
    $overlap_group_info->{'start'} = $segment->{'start'};
    $overlap_group_info->{'stop'} = $segment->{'stop'};
    $overlap_group_info->{'segments'} = \@overlap_array;
    push(@{$overlap_group_info->{'segments'}}, $segment);
}

#take care of last one in case of overlap extending to the end of the last chromosome
print_overlaps();

#if no segments, print none
if ($overlap_counter == 0) {
    print "No overlapping segments found\n";
}
