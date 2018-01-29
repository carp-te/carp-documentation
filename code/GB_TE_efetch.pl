##This script is used to retrieve sequences from NCBI based on following key words:
##"reverse transcriptase", "transposon", "repetitive element", "RNA-directed DNA polymerase",
##"pol protein", "non-LTR retrotransposon", "mobile element", "retroelement", "polyprotein", 
##"retrovirus", ("group-specific antigen" "gag"), "polymerase (pol)"

##Usage: perl GB_TE_efetch.pl fasta <output>; <output> is the file name used to store retrieved sequences.

#! /usr/bin/perl 
$db=$ARGV[0];   
@timeData = localtime(time);
print "This progarm starts at\t".join(' ', @timeData)."\n";
use Bio::DB::EUtilities;
# set optional history queue
my $factory = Bio::DB::EUtilities->new(-eutil => 'esearch',
                         -db => 'protein',
                         -term => '"reverse transcriptase" or "transposon" or "repetitive element" or "RNA-directed DNA polymerase" or "pol protein" or "non-LTR retrotransposon" or "mobile element" or "retroelement" or "polyprotein" or "retrovirus" or ("group-specific antigen" "gag") or "polymerase (pol)"',
                         -usehistory => 'y',
                         -rettype => $db);
 
my $count = $factory->get_count;
print "There are total $count\n items existing in protein databases.";
# get history from queue
my $hist = $factory->next_History || die 'No history data returned';
print "History returned\n";
# note db carries over from above
$factory->set_parameters(-eutil => 'efetch',
#                         -rettype => 'genbank',
                         -history => $hist);
my ($retmax, $retstart) = (500,0);
RETRIEVE_SEQS:
while ($retstart < $count) {
    $factory->set_parameters(-retmax => $retmax,
                            -retstart => $retstart);
    eval{
        $factory->get_Response(-file => "efetch.tmp");
    };
    if ($@) {
        die "Server error: $@.  Try again later" if $retry == 5;
        print STDERR "Server error, redo #$retry\n";
        $retry++ && redo RETRIEVE_SEQS;
    }
    open(TMP,"efetch.tmp");
    open(OUT,">>".$ARGV[1]);
    while(<TMP>){
	print OUT;
    }
    $retstart += $retmax;
}
close TMP; close OUT;
@timeData = localtime(time);
print "This progarm ends at\t".join(' ', @timeData)."\n"; 
