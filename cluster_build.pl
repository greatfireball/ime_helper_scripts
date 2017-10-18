#!/usr/bin/env perl

use strict;
use warnings;

use Graph::Undirected;

my $inputfile = shift @ARGV;
my $clustersize = shift @ARGV || 6000;

my %chromosomes = ();
my %pos_2_mirna = ();

open(FH, "<", $inputfile) || die "$!";
while (<FH>)
{
    chomp;
    s/"//g;

    my ($pos_number, $genome, $start, $stop, $chr, $strand, $mirna)=split(",", $_);

    my $mirna_pos_entry = {pos_number => $pos_number, genome => $genome, start => $start, stop => $stop, chr => $chr, strand => $strand, mirna => $mirna};
    push(@{$chromosomes{$chr}{$strand}}, $mirna_pos_entry);
    $pos_2_mirna{$pos_number} = $mirna_pos_entry;
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
	    my $g = Graph::Undirected->new;

	    for (my $i=0; $i<@{$chromosomes{$chr}{$strand}}-1; $i++)
	    {
		for (my $j=$i+1; $j<@{$chromosomes{$chr}{$strand}}; $j++)
		{
		    # add both nodes if they are not already present
		    foreach my $v ($i, $j)
		    {
			unless ($g->has_vertex($chromosomes{$chr}{$strand}[$v]{pos_number}))
			{
			    $g->add_vertex($chromosomes{$chr}{$strand}[$v]{pos_number});
			}
		    }

		    # use the smaller distance between start of one and stop of the other miRNA position
		    my $dist = abs($chromosomes{$chr}{$strand}[$i]{start}-$chromosomes{$chr}{$strand}[$j]{stop});
		    if (abs($chromosomes{$chr}{$strand}[$i]{stop}-$chromosomes{$chr}{$strand}[$j]{start}) < $dist)
		    {
			$dist = abs($chromosomes{$chr}{$strand}[$i]{stop}-$chromosomes{$chr}{$strand}[$j]{start});
		    }

		    # if the dist is smaller than the maximal allowed cluster distance, create a connection
		    if ($dist <= $clustersize)
		    {
			# add an edge between the nodes
			$g->add_edge($chromosomes{$chr}{$strand}[$i]{pos_number}, $chromosomes{$chr}{$strand}[$j]{pos_number});
		    }
		}
	    }

	    # get all connected components as clustes
	    my @clusters = $g->connected_components();

	    foreach my $cluster (@clusters)
	    {
		my @clustered_miRNAs = sort {$a->{start} <=> $b->{stop} || $a->{stop} <=> $b->{stop}} map {$pos_2_mirna{$_}} (@{$cluster});
		if (@clustered_miRNAs > 1)
		{
		    my $cluster_start = $clustered_miRNAs[0]{start};
		    my $cluster_stop  = $clustered_miRNAs[@clustered_miRNAs-1]{start};

		    printf "Cluster of miRNAs (border-size: %d) on chromosome %s(%s) (%d-%d) containing %d miRNAs: %s\n", $clustersize, $chr, $strand, $cluster_start, $cluster_stop, (@clustered_miRNAs+0), join(", ", map {$_->{mirna}} (@clustered_miRNAs));
		} else {
		    my $singleton = $clustered_miRNAs[0];
		    printf STDERR "Single miRNA cluster for miRNA %d on chromosome %s(%s) (%d-%d)\n", $singleton->{mirna}, $chr, $strand, $singleton->{start}, $singleton->{stop};
		}
	    }

	} else {
	    die "Should never happen";
	}
    }
}
