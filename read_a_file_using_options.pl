#!/usr/bin/perl

use strict;
use Getopt::Std;

our ($opt_i);
getopts('i:');
if (!$opt_i) {
    print STDERR "\nInput file name (-i) required\n\n\n";
}

open INFILE, "<", $opt_i or die "No such input file $opt_i";


my $first_line = 1; #used to skip the header;

#reads the file line by line

while (<INFILE>) {

    #skips the header
    if ($first_line) { 
	$first_line = 0;
	next;
    }

    # the variable $_ contains the current line.

    chomp $_; #strips the new line character from the current input

    # turns the current tab delimited line into an array of values
    my @row =  split('\t', $_);


    #here is where we actually process the line.



}

close INFILE;
