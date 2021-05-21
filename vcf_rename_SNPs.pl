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


    print $row[0]."\t".$row[1]."\t"."S".$row[0]."_".$row[1]."_".$row[4]."\t".$row[3]."\t".$row[4]."\t".$row[5];
    
    #ignore hets and anything not 0/0 or 1/1
    for (my $i=6; $i<=$#row;$i++ ) {
	my $cell = $row[$i];
	print "\t".$cell;
    }
    print "\n";
}

close INFILE; #close input file after all lines have been processed


