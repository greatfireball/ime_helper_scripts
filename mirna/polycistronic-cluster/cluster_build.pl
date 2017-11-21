#!/usr/bin/env perl

use strict;
use warnings;

use Graph::Undirected;

=pod

=head1 cluster_build.pl

A small program to generate miRNA clusters for our database based on a
table of miRNA positions.

=head1 SYNOPSIS

    # build a list of clusters with max distance 7500 nt from input.csv
    cluster_build.pl input.csv 7500

=head1 PARAMETER

=over 4

=item Input file

First parameter is the input csv file.

=item cluster size

Second parameter is the maximum distance between two miRNAs still
belonging to the same cluster. If no parameter is provided, 6000 is
assumend. If multiple clustersizes should be testes, just put them
comma/space separated as second parameter.

=back

=head1 INPUT

A csv file with comma as field separator containing the following
information:

=over 4

=item position_entry_id [mandatory, unique]

=item genome_id

=item start_position [mandatory]

=item stop_position [mandatory]

=item chromosome [mandatory]

=item strand [mandatory]

=item miRNA_precursor_ID [mandatory]

=back

Quotes will be stripped out.

=head1 OUTPUT

Output will be printed to STDOUT and additional information are
printed to STDERR. The csv content of STDOUT contains the following
columns:

=over 4

=item cluster_id

=item position_entry_id

=item maximal_distance

=back

The columns are comma separated. position_entry_id will be returned
due to a single precursor might have multiple copies inside the
genome.

=head1 AUTHOR

This script was written by Frank FE<ouml>rster C<frank.foerster@ime.fraunhofer.de>.

=head1 CHANGELOG

=over 4

=item 2017-10-18 v0.1.2

First working version.

=item 2017-10-18 v0.1.3

Updated documentation/allows multiple cluster sizes

=back

=head1 COPYRIGHT AND LICENCE

MIT License

Copyright (c) 2017 Frank FE<ouml>rster

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

=cut

use version 0.77; our $VERSION = version->declare("v0.1.3");

unless (@ARGV)
{
    die "You need to specify at least an input file. For further information see perldock cluster_build.pl.";
}

my $inputfile = shift @ARGV;
my @clustersizes = (@ARGV);

unless (@clustersizes)
{
    @clustersizes = (6000);
}

@clustersizes = split(/,/, join(",", @clustersizes));

my %chromosomes = ();
my %pos_2_mirna = ();
my $clusterid = 0;

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


foreach my $clustersize (@clustersizes)
{

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

			printf STDERR "Cluster of miRNAs (border-size: %d) on chromosome %s(%s) (%d-%d) containing %d miRNAs: %s\n", $clustersize, $chr, $strand, $cluster_start, $cluster_stop, (@clustered_miRNAs+0), join(", ", map {$_->{mirna}} (@clustered_miRNAs));

			# print the output for the cluster table
			$clusterid++;
			foreach my $position_id (@{$cluster})
			{
			    print join(",", ($clusterid, $position_id, $clustersize)), "\n";
			}
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
}
