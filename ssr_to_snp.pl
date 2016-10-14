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

my $first_line = 1; #used to skip the header;

#reads the file line by line
while (<INFILE>) {

    # the variable $_ contains the current line.
    chomp $_; #strips the new line character from the current input

    #skips the header
    if ($first_line) { 
	#also print the header to the output
	print "$_\n";
	$first_line = 0;
	next;
    }

    #turns the current tab delimited line into an array of values
    my @row =  split('\t', $_);

    ###################################
    #here is where we actually process the row
    #the objective in this case is to first get allele frequencies
    ###################################

    #read and removes the first "cell" of the row
    #the script could be modified to store and reprint additional columns
    my $marker_name = shift @row;

    #todo: check for presence of a slash and non-numeric data

    my @all_alleles; #will be a list of alleles used to calculate frequency

    #split each cell by slash
    foreach my $cell (@row) {
	my @split_cell = split('\/', $cell);
	push @all_alleles, @split_cell; #adds the current array of alleles to the total list
    }

    #creates a hash to store alleles and counts
    my %allele_counts;

    #count all of the alleles for this marker
    foreach my $allele (@all_alleles) {

	#skip non-integers
	unless ($allele =~ /^\d+?$/) {
	    next;
	}

	if ( $allele_counts{$allele} ) { #checks if this allele has been seen yet
	    $allele_counts{$allele}++; #adds one to the allele count
	} else {
	    $allele_counts{$allele} = 1; #sets the allele count if first time seen
	}
    }

    #make a list of alleles sorted by order NOTE: this sort is only for numberic values
    #sorted by values (frequency) not keys
    my @allele_frequency_sort = sort { $allele_counts{$a} <=> $allele_counts{$b} } keys %allele_counts;
    

    ####################################################
    # now that we have the alleles sorted by frequency we will assign fake nuckeotides
    ####################################################

    #make an hash to reassign bases to numeric allele calls
    my %allele_decoder;

    #make an array of fake nucleotide codes
    my @nucleotides = qw( A C G T + -);  #must be in order we wish to them from max to min


    #for each fake nucleotide this finds the next most frequent allele size
    foreach my $base (@nucleotides) {

	#check for empty list in case all alleles have already been assigned a fake nucleotide
	if (scalar @allele_frequency_sort < 1) {
	    next; #just leave the foreach loop if this happens.  there is nothing left to do
	}

	#gets the most frequent allele size from the sorted list and removes it from the list
	#so that the next fake nucleotide will be assigned to next most frequent size
	my $current_allele = shift @allele_frequency_sort;

	#warn when two alleles have the same frequency.
	#this won't hurt anything but the assignment my be different each time the script is run
	if (scalar @allele_frequency_sort > 1 && $allele_counts{$current_allele} == $allele_counts{$allele_frequency_sort[0]}) {
	    print stderr "warning $marker_name allele $current_allele and allele ".$allele_frequency_sort[0]." have equal frequency\n";
	}

	#set the decoder so that the current allele size is assigned the fake nucleotide
	$allele_decoder{$current_allele} = $base;

	#just helps the user see what is happening.
	#this could be sent to an optional file if recording the allele size to base coding
	print stderr "$marker_name allele ".$current_allele." has been set to nucleotide $base\n";
      }

    #set all remaining alleles to the "-" code
    #this could also be sent to a file somewhere
    foreach my $remaining_allele (@allele_frequency_sort) {
	$allele_decoder{$remaining_allele} = "-";
	print stderr "$marker_name allele ".$remaining_allele." has been collapsed and set to nucleotide -\n";
    }


    ####################################################
    # now that we have the decoder we can reprint the row
    ####################################################

    #marker name is unchanged
    #the script could be modified to store and reprint additional columns
    print "$marker_name";

    #go cell by cell through the line to reprint
    foreach my $cell (@row) {

	print "\t"; #tab delimited

	my @split_cell = split('\/', $cell); # split on / character

	#go through / delimited alleles for the individual cell
	#note this also should work with polyploids but has only been tested with diploids
	foreach my $old_allele (@split_cell) {

	    #lookup current allele in decoder and print its corresponding fake nucleotide
	    if ($allele_decoder{$old_allele}) {
		print $allele_decoder{$old_allele};
	    }
	    #if the current allele is not in the decoder print an N
	    #this should take care of various forms of missing data coded as ? . or just empty
	    #this behavior could be changed to pass through the non-numeric alleles as is
	    else {
		#print $old_allele; #uncomment to pass through non-numeric and comment out below
		print "N"; 
	    }
	}
    }
    print "\n"; #ends the line
}

close INFILE; #close input file after all lines have been processed


