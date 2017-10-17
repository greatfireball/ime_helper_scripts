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

    push(@{$chromosomes{$chr}{$strand}}, {start => $start, stop => $stop, chr => $chr, strand => $strand, mirna => $mirna});
}
close(FH) || die "$!";

foreach my $chr (keys %chromosomes)
{
    my $num_mirnas_plus = (exists $chromosomes{$chr}{'+'}) ? @{$chromosomes{$chr}{'+'}} : 0;
    my $num_mirnas_minus = (exists $chromosomes{$chr}{'-'}) ? @{$chromosomes{$chr}{'-'}} : 0;
    if ($num_mirnas_plus <= 1)
    {
	print STDERR "Skipping chromosome $chr plus strand, due to its single miRNA content.\n";
    }
    if ($num_mirnas_minus <= 1)
    {
	print STDERR "Skipping chromosome $chr minus strand, due to its single miRNA content.\n";
    }

}
