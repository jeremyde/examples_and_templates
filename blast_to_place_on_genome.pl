#!/usr/bin/perl

=pod

=head1 NAME

blast_to_place_on_genome.pl

=head1 SYNOPSIS

blast_to_place_on_genome.pl -d database.fa -f sequences.fa -l 500 -p outprefix 

=head1 DESCRIPTION

This script locates the genomic region that is most likely matching a sequence such as an RFLP probe


=head1 LICENSE

  Same as Perl.

=head1 AUTHORS

  Jeremy D. Edwards <jde22@cornell.edu>

=cut

use strict;
use warnings;
use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::DB::Fasta;
use Bio::Tools::Run::StandAloneBlast;

our ($opt_d, $opt_f, $opt_l, $opt_p, $opt_h);
getopts('d:f:l:p:h');
if ($opt_h){
    help();
    exit;
}
if ((!$opt_d) || (!$opt_f) || (!$opt_l) || (!$opt_p)) {
    print STDERR "\nMissing required options\n\n\n";
    help();
}

# User options
my $database = $opt_d;    # FASTA file formatted for BLAST with formatdb
my $mfasta_in = $opt_f;      # Multi FASTA file of sequences to place
my $max_length = $opt_l;  # maximum distance between HSPs in bp
my $out_found = $opt_p.".found.csv";       #output file for matched sequences
my $out_notfound = $opt_p.".notfound.csv"; #output file for unmatched sequences

# built in
my $evalue=0.00001;

# Data storage
my @found_sequences;
my @notfound_sequences;

# open output files
open(FOUND,"> $out_found")||die("can't open $out_found: $!\n");
open(NOTFOUND,"> $out_notfound")||die("can't open $out_notfound: $!\n");

# open multi FASTA input file
open(INFILE,$mfasta_in)||die("can't open $mfasta_in: $!\n");
close(INFILE);

my $seq_in;
eval {
   $seq_in   = Bio::SeqIO->new(
                               -format => 'fasta',
                               -file   => $mfasta_in,
                               );
}; 

if( $@ ) {
   print "Was not able to open file $mfasta_in";
   print "Full error is $@ ";
   exit(-1);
}

# print column headers in output files
print FOUND "ID\tChrom\tStart\tStop\tStrand\tevalue\tHSPs\n";
print NOTFOUND "ID\n";


    # setup BLAST
    my $factory = Bio::Tools::Run::StandAloneBlast->new('program'  => 'blastn',
							'database' => $database,
							'F' => 'F',
							'W' => 20,
							'a' => 2,
							_READMETHOD => "Blast",
							'e' => $evalue
	);


while (my $seq = $seq_in->next_seq()) {

    # Run BLAST
    my $blast_report = $factory->blastall($seq);
    my $result = $blast_report->next_result; #only one BLAST search so only one result
    my $hit = $result->next_hit(); #only using best hit
    my $subject_start;
    my $subject_end;
    my $hsp_found;
    my $best_hsp_evalue;
    my $hsp_count = 1;
    unless ($hit) { #deal with no hits
	print NOTFOUND $seq->id()."\n";
	next;
    }
    if ($hit->significance() > $evalue) { #deal with nonsignificant hits
	print NOTFOUND $seq->id()."\n";
	next;
    }
    while (my $hsp = $hit->next_hsp()) {
	if ($hsp->evalue() <= $evalue) {

	    print "id: ".$seq->id()."\n";
	    print "eval: ".$hit->significance()."\n";
	    print "h_eval: ".$hsp->evalue()."\n";
	    print "chr: ".$hit->name()."\n";
	    print "strand: ".$hit->strand('subject')."\n";
	    print "h_strand: ".$hsp->strand('subject')."\n";
	    print "start: ".$hsp->start('subject')."\t";
	    print "end: ".$hsp->end('subject')."\n\n";
	    unless ($hsp_found) { 
		$hsp_found = 1;
		$subject_start = $hsp->start('subject');
		$subject_end = $hsp->end('subject');
		$best_hsp_evalue = $hsp->evalue();
	    }		
	    
	    if ($hsp->start('subject') < $subject_start) {
		if ($subject_start - $hsp->start('subject')  <= $max_length) {
		    $subject_start = $hsp->start('subject');
		    $hsp_count++;
		} else {
		    next;
		}
	    }
	    if ($hsp->end('subject') > $subject_end) {
		if ($hsp->end('subject') - $subject_end  <= $max_length) {
		    $subject_end = $hsp->end('subject');
		    $hsp_count++;
		} else {
		    next;
		}
	    }
	} 
    }
    if ($hsp_found) {
	print FOUND $seq->id()."\t".$hit->name()."\t".$subject_start."\t".$subject_end."\t".$hit->strand('subject')."\t".$best_hsp_evalue."\t".$hsp_count."\n";
    } else {
	print NOTFOUND $seq->id()."\n";
    }
}
$factory->cleanup();



# my $linenum=0;
# while (<INFILE>) {
#     chomp $_;
#     if ($linenum == 0) { # skip the header
# 	$linenum++;
# 	next;
#     }
#     my $temp = $_;
#     $temp =~ s/"//g;    # strip the quotes
#     $temp =~ s///g;   # strip those weird ^M characters that show up as line terminators

#     #This is the order in which columns are expected in the input file:
#     my @line = split(/\t/,$temp);
#     my $marker_name = $line[$ID_col];
#     my $primer1 = $line[$primer1_col];
#     my $primer2 = $line[$primer2_col];

#     # get rid of whitespace in any of these
#     if (($marker_name ne "") && ($primer1 ne "") && ($primer2 ne "")) {
# 	$marker_name =~ s/\s//g;
# 	$primer1 =~ s/\s//g;
# 	$primer2 =~ s/\s//g;
#     }

#     # skip row if any of the fields are empty or invalid
#     if (!(($marker_name ne "") && ($primer1 ne "") && ($primer2 ne "") && ($primer1 =~ /\A[acgtACGT]*\z/i) && ($primer1 =~ /\A[acgtACGT]*\z/i))) {
# 	print "Skipping amplicon $marker_name line $linenum\n";
# 	next;
#     }

#     # store concatenated string of primers
#     my $string = $primer1.$spacer.$primer2;
#     $markers{$marker_name}=$string;

#     my $i = 0; # increment counter for evalues
#     while ((not exists $positive_match{$marker_name}) && ($i <= $#evalues)) {

# 	# create sequence object from primer string
# 	my $query = Bio::Seq->new(-display_id => $marker_name, -seq => $string, -alphabet => 'dna');

# 	# setup BLAST
# 	my $factory = Bio::Tools::Run::StandAloneBlast->new('program'  => 'blastn',
# 							    'database' => $database,
# 							    'F' => 'F',
# 							    _READMETHOD => "Blast",
# 							    'e' => $evalues[$i]
# 	    );
# 	# Run BLAST
# 	my $blast_report = $factory->blastall($query);
# 	my $result = $blast_report->next_result; #only one BLAST search so only one result

# 	my $strandsum = 0;
# 	my @hsp_locations = ();
# 	my $match_flag = 0; # is this used?

# 	my %forwardhash;	# HSPs on strand 1
# 	my %reversehash;	# HSPs on strand -1

# 	while (my $hit = $result->next_hit()) { # loop through all hits

# 	    print REPORT "$evalues[$i]\t".$query->display_id()."\t",
# 	    $hit->accession,"\t",
# 	    $hit->num_hsps,"\t";
# 	    undef @hsp_locations;

# 	    $strandsum = 0;
# 	    foreach my $item ($hit->hsps) {
# 		if ($item->strand('hit') eq 1) {
# 		    $forwardhash{$item->start('hit')}=$item;
# 		} elsif ($item->strand('hit') eq -1) {
# 		    $reversehash{$item->end('hit')}=$item;
# 		}

# 		$strandsum = $strandsum + $item->strand('hit');
# 		push (@hsp_locations,$item->start('hit'),$item->end('hit'));
# 		print REPORT "\t",$item->start('query'),"-",$item->end('query'),
# 		"\t",$item->strand('hit'),
# 		"\t",$item->evalue,
# 		"\t",$item->start('hit'),"-",$item->end('hit'),"\t",
# 	    }		# while parsing hsps
# 	    print REPORT "\n";

# 	    foreach my $fwd (sort {$a<=>$b} keys %forwardhash) {
# 		foreach my $rev (sort {$b<=>$a} keys %reversehash) {
# 		    if (($fwd <= $rev) && ($rev-$fwd <= $max_length)) {

# 			$num_matches{$query->display_id()}++;
# 			my $IDandnum_matches = $query->display_id()."\t".$num_matches{$query->display_id()};
# 			my $size = $rev-$fwd;
# 			$matches{$IDandnum_matches} =
# 			    $query->display_id()."\t".		   # Marker
# 			    $hit->accession."\t".		   # Chrom
# 			    $size."\t".			   # size
# 			    "$fwd"."-"."$rev"."\t".		   # region
# 			    $forwardhash{$fwd}->evalue."\t". # left primer eval
# 			    $forwardhash{$fwd}->start('query')."\t". # left primer start position
# 			    $forwardhash{$fwd}->end('query')."\t". # left primer size
# 			    $reversehash{$rev}->evalue."\t". # right primer eval
# 			    $reversehash{$rev}->start('query')."\t". # right primer start position
# 			    $reversehash{$rev}->end('query'). # right primer size
# 			    "\t$string"; # query string

# 			$positive_match{$query->display_id()}=$hit->accession;
# 			$match_flag = 1;  # is this used?
# 		    }
# 		}
# 	    }

# 	    undef %forwardhash;
# 	    undef %reversehash;

# 	}			# while parsing hits
# 	$i++;
#     }		 # while incrementing $i until we get a positive match
#     $linenum++;
# }
# close (INDEX);
# close (REPORT);

# foreach my $marker (sort keys %num_matches) {
#   if ($num_matches{$marker} <= $max_matches) {

#     foreach my $match (sort keys %matches) {
#       if ($match =~ /$marker\t/) {
# 	print MATCHES "$num_matches{$marker}\t$matches{$match}\n";
#       }
#     }
#   } else {
#     foreach my $match (sort keys %matches) {
#       if ($match =~ /$marker\t/) {
# 	print UNMATCHED "$num_matches{$marker}\t$matches{$match}\n";
#       }
#     }
#   }
# }

# foreach my $marker (sort keys %markers) {
#   if (not exists $positive_match{$marker}) {
#     print  UNMATCHED "0\t$marker\t\t\t\t\t\t\t\t\t\t$markers{$marker}\n";
#   }
# }
# close(MATCHES);
# close(UNMATCHED);




sub help {
  print STDERR <<EOF;
  $0:

    Description:

      This script locates the genomic region that is most likely matching a sequence such as an RFLP probe

    Usage:

      blast_to_place_on_genome.pl -d database.fa -f sequences.fa -l 500 -p outprefix 

    Flags:

      -d <database_file>                 FASTA database file. Must be formatted with formatdb. (mandatory)

      -f <fasta_file>                    FASTA file of sequences to place

      -l <max_length>                    Maximum distance between HSPs in bp (mandatory)

      -p <output_file_prefix>            Output file prefix (mandatory)

      -h <help>                          Help

EOF
exit (1);
}
