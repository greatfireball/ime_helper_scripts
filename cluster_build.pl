#!/usr/bin/env perl

use strict;
use warnings;

my $inputfile = shift @ARGV;

my %chromosomes = ();

open(FH, "<", $inputfile) || die "$!";
while (<FH>)
{
    chomp;
    s/"//g; 

    my (undef, undef, $start, $stop, $chr, $strand, $mirna)=split(",", $_); 

    push(@{$chromosomes{$chr}}, {start => $start, stop => $stop, chr => $chr, strand => $strand, mirna => $mirna});
}
close(FH) || die "$!";
