#!/usr/bin/perl

use strict;
use Getopt::Std;

#reading options
our ($opt_i, $opt_h);
getopts('i:h:');
if (!$opt_i) {
    print STDERR "\nInput file name (-i) required\n\n\n";
}

my $het_limit = $opt_h;

#open the file or die
open INFILE, "<", $opt_i or die "No such input file $opt_i";


#reads the file line by line
while (<INFILE>) {


    # the variable $_ contains the current line.
    chomp $_; #strips the new line character from the current input

    if ($_ =~ /^#/) {
	print "$_\n";
	next;
    }

    #turns the current delimited line into an array of values
    my @row =  split(' ', $_);

    my $het_count = 0;
    my $total_count = 0;

    foreach my $cell (@row) {
	if ($cell =~ /\//) {
	    $total_count++;
	    if ($cell eq '0/1') {
		$het_count++;
	    }
	}
    }

    if ($total_count > 0) {
	my $het_ratio = $het_count/$total_count;
	if ($het_ratio < $het_limit) {
	    print "$_\n";
	}
    }
}

close INFILE; #close input file after all lines have been processed


