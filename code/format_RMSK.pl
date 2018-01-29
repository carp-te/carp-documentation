#!/bin/usr/perl -w
use strict;

=head1 SYNOPSIS
format_RMSK.pl <RMSKfile> > <output>
=cut

my ($filename) = @ARGV;
open(IN,"$filename");

# Change RepeatMasker output into CENSOR map similar results
while(<IN>){
    chomp;
    my @data = split(" ",$_);
    if($data[8] eq '+'){
	if ($data[11] =~ /\([\w]+\)/){
	    print "$data[4]\t$data[5]\t$data[6]\t$data[9]\t$data[13]\t$data[12]\td\t0\t0\t$data[0]\n";
	}
	else {
	    print "$data[4]\t$data[5]\t$data[6]\t$data[9]\t$data[11]\t$data[12]\td\t0\t0\t$data[0]\n";
	}
	
    }
    if($data[8] eq 'C'){
	if ($data[11] =~ /\([\w]+\)/){
	    print "$data[4]\t$data[5]\t$data[6]\t$data[9]\t$data[13]\t$data[12]\tc\t0\t0\t$data[0]\n";
	}
	else {
	    print "$data[4]\t$data[5]\t$data[6]\t$data[9]\t$data[11]\t$data[12]\tc\t0\t0\t$data[0]\n";
	}
    }

}
