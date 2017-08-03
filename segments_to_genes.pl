#!/usr/bin/perl

#inputvariables segment_file gene_file(sorted)

use strict;
use Getopt::Std;
use DateTime qw( );

#reading options
our ($opt_s, $opt_g);
getopts('s:g:');
if (!$opt_s) {
    print STDERR "\nSegment input file name (-i) required\n\n\n";
}
if (!$opt_g) {
    print STDERR "\nGene list file name (-i) required\n\n\n";
}

open GENEFILE, "<", $opt_g or die "No such input file $opt_g";
open SEGFILE, "<", $opt_s or die "No such input file $opt_s";

my @descriptions;
my @chromosomes;
my @starts;
my @ends;
my @strands;

while (<GENEFILE>) {
    chomp $_; #strips the new line character from the current input
    my @row =  split('\t', $_);
    push @chromosomes, $row[0];
    push @starts, $row[1];
    push @ends, $row[2];
    push @strands, $row[3];
    push @descriptions, $row[4];
}
close GENEFILE;

my @segmentchromosomes;
my @segmentstarts;
my @segmentends;
my @segmentpeaks;
my @segmentpeakvals;

my $firstline = 1;
while (<SEGFILE>) {
    #skip header
    if ($firstline == 1) {
	$firstline = 0;
	next;
    }
    chomp $_; #strips the new line character from the current input
    my @row =  split('\t', $_);
    push @segmentchromosomes, $row[0];
    push @segmentstarts, $row[1];
    push @segmentends, $row[2];
    push @segmentpeaks, $row[3];
    push @segmentpeakvals, $row[4];
}
close SEGFILE;


print "Genes found within chromosome segments\n";
print "Report for file $opt_s using genes found in file $opt_g\n";
my $dt = DateTime->now( time_zone => 'local' );
print "Date: ".$dt->strftime("%m/%d/%Y")."\n";

for (my $i = 0; $i < scalar(@segmentchromosomes); $i++) {
    my @genelist;
    for (my $j = 0; $j < scalar(@chromosomes); $j++) {
	if ($chromosomes[$j] eq $segmentchromosomes[$i]) {
	    my $genestart;
	    my $geneend;
	    # get orientation.  could use strand instead
	    if ($starts[$j] < $ends[$j]) {
		$genestart = $starts[$j];
		$geneend = $ends[$j];
	    } else {
		$genestart = $ends[$j];
		$geneend = $starts[$j];
	    }
	    # stop searching once past the gene.  List assumed sorted
	    if ($genestart > $segmentends[$i]) {
		last;
	    }
	    if ($geneend < $segmentstarts[$i]) {
		next;
	    }
	    if ((($genestart > $segmentstarts[$i]) && ($genestart < $segmentends[$i])) || (($geneend > $segmentstarts[$i]) && ($geneend < $segmentends[$i]))) {
		my %gene;
		$gene{'chr'} = $chromosomes[$j];
		$gene{'start'} = $starts[$j];
		$gene{'end'} = $ends[$j];
		$gene{'strand'} = $strands[$j];
		$gene{'desc'} = $descriptions[$j];
		if ($genestart < $segmentpeaks[$i] && $geneend > $segmentpeaks[$i]) {
		    $gene{'dist'} = 0;
		} elsif ($geneend < $segmentpeaks[$i]) {
		    $gene{'dist'} = $geneend - $segmentpeaks[$i];
		} else {
		    $gene{'dist'} = $genestart - $segmentpeaks[$i];
		}
		push @genelist, \%gene;
	    }
	}
	       }
    my @peak_in_gene;
    my @before_peak;
    my @after_peak;
    #make three lists of genes sorted by distance from peak snp
    foreach my $gene (@genelist) {
	if ($gene->{'dist'} == 0) {
	    push @peak_in_gene, $gene;
	}
	if ($gene->{'dist'} > 0) {
	    push @after_peak, $gene; #add to end
	}
	if ($gene->{'dist'} < 0) {
	    unshift @before_peak, $gene; #add to front
	}
    }
    sub commify {
	my $text = reverse $_[0];
	$text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
	return scalar reverse $text
    }
    sub PrintGenes {
	my @printlist = @_;
	print "distance\tchr\tstart\tend\tstrand\tdescription\n";
	foreach my $gene (@printlist) {
	    print commify($gene->{'dist'})."\t";
	    print $gene->{'chr'}."\t";
	    print commify($gene->{'start'})."\t";
	    print commify($gene->{'end'})."\t";
	    print $gene->{'strand'}."\t";
	    print $gene->{'desc'}."\n";
	}
	print "\n";
    }
    print "\n**********************************************************************************************************************************\n";
    print "\nSegment:\n";
    print "$segmentchromosomes[$i]\t".commify($segmentstarts[$i])."\t".commify($segmentends[$i])."\n";
    print "Peak SNP at position ".commify($segmentpeaks[$i])." p-value ".commify($segmentpeakvals[$i])."\n\n";
    if (scalar(@peak_in_gene) > 0) {
	print "Genes spanning peak SNP:\n";
	PrintGenes(@peak_in_gene);
    }
    if (scalar(@before_peak) > 0) {
	print "Genes before peak SNP:\n";
	PrintGenes(@before_peak);
    }
    if (scalar(@after_peak) > 0) {
	print "Genes after peak SNP:\n";
	PrintGenes(@after_peak);
    }
}


