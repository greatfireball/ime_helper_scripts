#!/usr/bin/env perl

use strict;
use warnings;

my %gff = ();
my $counter = 0;

my $file = "manduca_genomic_regions.gff";

open(FH, "<", $file) || die "$!\n";

while (<FH>)
{
    chomp();
    next if (/^#/);

    my ($chr, $source, $type, $start, $stop, $score, $strand, $phase, $attribute) = split(/\t/, $_);

    # split the attributes
    my %attributes = ();
    foreach my $attr (split(/;/, $attribute))
    {
	my ($key, $value) = split(/=/, $attr);
	if (defined $key && defined $value)
	{
	    push(@{$attributes{$key}}, split(/,/, $value));
	} else {
	    warn "Error on parsing the attribute part '$attr'\n";
	}
    }

    my $name = sprintf("generated%05d", $counter++);
    if (exists $attributes{ID})
    {
	if (@{$attributes{ID}}==1)
	{
	    $name = $attributes{ID}[0];
	} else {
	    warn "Multiple entries for ID in line '$_'\n";
	}
    }

    push(@{$gff{$chr}}, [$source, $type, $start, $stop, $score, $strand, $phase, \%attributes]);
}

close(FH) || die "$!\n";

warn "Sorting information\n";

foreach my $chr (keys %gff)
{
    @{$gff{$chr}} = sort {$a->[2] <=> $b->[2] || $a->[3] <=> $b->[3] || $a->[1] cmp $b->[1]} @{$gff{$chr}};
}

warn "Finished\n";
