#!/usr/bin/perl

use strict;
use Getopt::Std;

#reading options
our ($opt_i);
getopts('i:');
if (!$opt_i) {
    print STDERR "\nInput file name (-i) required\n\n\n";
}


#open the file or die
open INFILE, "<", $opt_i or die "No such input file $opt_i";

my @acclist;
my $firstline = 1;


#reads the file line by line
while (<INFILE>) {


    # the variable $_ contains the current line.
    chomp $_; #strips the new line character from the current input

    if ($_ =~ /^##/) {
	print "$_\n";
	next;
    }

    
    #turns the current delimited line into an array of values
    my @row =  split(' ', $_);

    if ($firstline) {
	for (my $i=9; $i<=$#row;$i++ ) {
	    my @accname = split('_',$row[$i]);
	    push @acclist, $accname[1];
	}
	print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
	foreach my $acc (@acclist) {
	    print "\t$acc";
	}
	print "\n";
	$firstline = 0;
	next;
    }	

    print "$_\n";
}

close INFILE; #close input file after all lines have been processed


