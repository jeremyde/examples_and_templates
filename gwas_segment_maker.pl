#!/usr/bin/perl

#inputvariables file sig_col_name sig_thresh extend_sig_thresh distance join_distance
#inputfile SNP	CHR	BP     	P	qvalue
#outputfile CHR BP Start Stop Peak Peak_val

use strict;
use Getopt::Std;
use Math::Complex;
use Scalar::Util qw(looks_like_number);

#reading options
our ($opt_i, $opt_c, $opt_t, $opt_e, $opt_d, $opt_j);
getopts('i:c:t:e:d:j:');
if (!$opt_i) {
    print STDERR "\nInput file name (-i) required\n\n\n";
}
#todo: check other opts

#open the file or die
open INFILE, "<", $opt_i or die "No such input file $opt_i";

my $first_line = 1; #used to recognize the header;
my @header_row;
my $chr_col;
my $bp_col;
my $sig_col;
my $peak_snp;
my $peak_snp_value;

my @segments_discovered;

my $segment_open = 0;
my $segment_start;
my $segment_end;
my $segment_last;
my $current_chromosome;

#reads the file line by line
while (<INFILE>) {

    # the variable $_ contains the current line.
    chomp $_; #strips the new line character from the current input

    #skips the header
    if ($first_line) { 
	@header_row =  split('\t', $_);

	my $col_count = 0;
	foreach my $col_name (@header_row) {
	    if ($col_name eq "Chr") {
		$chr_col = $col_count;
	    }
	    if ($col_name eq "Pos") {
		$bp_col = $col_count;
	    }
	    if ($col_name eq $opt_c) {
		$sig_col = $col_count;
	    }
	    $col_count++;
	}
	#todo: check that all required columns are found or die

	$first_line = 0;
	next;
    }

    my @row =  split('\t', $_);

    #todo: add check for new chromosome


    #skip missing data
    if ($row[$sig_col] eq 'NaN') {
	next;
    }

    #check for new chromosome
    if ($current_chromosome) {
	if ($row[$chr_col] ne $current_chromosome) {
	    print stderr "Finished chr $current_chromosome, starting chr ".$row[$chr_col]."\n";
	    $current_chromosome = $row[$chr_col];
	    if ($segment_open == 1) {
		$segment_end = $segment_last;
		my %segment_hash;
		$segment_hash{'chr'} = $row[$chr_col];
		$segment_hash{'start'} = $segment_start;
		$segment_hash{'end'} = $segment_end;
		$segment_hash{'peak'} = $peak_snp;
		$segment_hash{'peak_val'} = $peak_snp_value;
		push @segments_discovered, \%segment_hash;
		$segment_open = 0;
	    }
	}
    }

    # open a new segment for significant SNPs if not already open
    if ($segment_open == 0) {
	#skip non significant rows
	if (-log10($row[$sig_col]) > $opt_t) {
	    $segment_open = 1;
	    $peak_snp = $row[$bp_col];
	    $peak_snp_value = $row[$sig_col];
	    $segment_start = $row[$bp_col];
	    $segment_last = $row[$bp_col];
	    $current_chromosome = $row[$chr_col];
	}
	next;
    }
    #segment is open
    if ($row[$bp_col] > $segment_last + $opt_d) {
	#close segment
	$segment_end = $segment_last;
	my %segment_hash;
	$segment_hash{'chr'} = $row[$chr_col];
	$segment_hash{'start'} = $segment_start;
	$segment_hash{'end'} = $segment_end;
	$segment_hash{'peak'} = $peak_snp;
	$segment_hash{'peak_val'} = $peak_snp_value;
	$segment_hash{'merged'} = "single";
	push @segments_discovered, \%segment_hash;
	$segment_open = 0;
	print stderr "Segment saved\n";
	next;
    }
    #segment is open and not past extend distance
    if (-log10($row[$sig_col]) > $opt_e) {
	$segment_last = $row[$bp_col];
	if (-log10($row[$sig_col]) > $peak_snp_value) {
	    $peak_snp = $row[$bp_col];
	    $peak_snp_value = $row[$sig_col];
	}
	next;
    }
}

close INFILE;

my %scaffold_segments;
my $first_segment = 1;
my $last_segment;
my $scaf_chr;
my $scaf_start;
my $scaf_end;
my $scaf_peak_snp;
my $scaf_peak_snp_val;
my $last_segment_printed = 0;

#print the header to the output
print "CHR\tStart\tStop\tPeak\tPeak_val\tMerged\n";

sub write_scaffold {
    #print $last_segment->{'chr'}."\t".$last_segment->{'start'} - $opt_d."\t".$last_segment->{'end'} + $opt_d."\t".$last_segment->{'peak'}."\t".$last_segment->{'peak_val'}."\t".$last_segment->{'merged'}."\n";
    my $adj_start = $last_segment->{'start'} - $opt_d;
    my $adj_end = $last_segment->{'end'} + $opt_d;
    print $last_segment->{'chr'}."\t".$adj_start."\t".$adj_end."\t".$last_segment->{'peak'}."\t".$last_segment->{'peak_val'}."\t".$last_segment->{'merged'}."\n";
}

foreach my $segment (@segments_discovered) {
    if ($first_segment == 1) {
	$current_chromosome = $segment->{'chr'};
	$first_segment = 0;
	$last_segment = $segment;
	next;
    }
    $last_segment_printed = 1;

    #catch and print joined segments at end of chromosome
    if ($segment->{'chr'} ne $last_segment->{'chr'}) {
	if ($last_segment_printed == 0) {
	    write_scaffold();
	    $last_segment_printed = 1;
	}
    }

    #see if within join distance
    if (($segment->{'chr'} eq $last_segment->{'chr'}) && ($last_segment->{'end'} + $opt_j > $segment->{'start'})) {
	print stderr "Merging adjacent scaffolds :".$last_segment->{'chr'}." ".$last_segment->{'start'}." ".$last_segment->{'end'}."\n";
	$last_segment->{'end'} = $segment->{'end'};
	#$last_segment->{'end'} = 0;
	$last_segment->{'merged'} = "merged";
	if ($last_segment->{'peak_val'} > $segment->{'peak_val'}) {
	    $last_segment->{'peak_val'} = $segment->{'peak_val'};
	    $last_segment->{'peak'} eq $segment->{'peak'};
	}
	print stderr "New end is :".$last_segment->{'end'}."\n";
	$last_segment_printed = 0;
	#write_scaffold();
	next;
 
    } else {
	#if there is a previous merged segment being build, print it now
	#if ($last_segment->{'merged'} eq "merged") {
	#    write_scaffold();
	#}
	#print the current segment
	write_scaffold();
	$last_segment = $segment;
	$last_segment_printed = 1;
    }
}

#write last scaffold
write_scaffold();

# todo: collect and print stats
    


