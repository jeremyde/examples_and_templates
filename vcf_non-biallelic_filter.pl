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

my $previous_line ;
my $previous_chromosome;
my $previous_position;
my $multi_on = 0;
my $first_line = 1;

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

    if ($first_line == 1) {
	$first_line = 0;
	$previous_line = $_;
	$previous_chromosome = $row[0];
	$previous_position = $row[1];
	next;
    }


    if (($row[0] eq $previous_chromosome ) && ($row[1] eq $previous_position)) {
	$multi_on = 1;
	next;
    }

    if ($multi_on == 0) {
	print "$previous_line\n";
    }

    $previous_line = $_;
    $previous_chromosome = $row[0];
    $previous_position = $row[1];
    $multi_on = 0;
    
}

if ($multi_on == 0) {
    print "$previous_line\n";
}

close INFILE; #close input file after all lines have been processed


