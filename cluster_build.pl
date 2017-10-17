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

foreach my $chr (sort keys %chromosomes)
{
    foreach my $strand (qw(+ -))
    {
	if (! exists $chromosomes{$chr}{$strand})
	{
	    print STDERR "Skipping chromosome $chr($strand), due to it contains no miRNA.\n";
	}
	elsif (@{$chromosomes{$chr}{$strand}} <= 1)
	{
	    print STDERR "Skipping chromosome $chr($strand), due to its single miRNA content.\n";
	} elsif (@{$chromosomes{$chr}{$strand}} > 1)
	{
	} else {
	    die "Should never happen";
	}
    }
}
