#!/usr/bin/env perl

use strict;
use warnings;

my $inputfile = shift @ARGV;
my $clustersize = shift @ARGV || 6000;

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
	    my @dna_strand = ();

	    # mark cluster for each mirna
	    foreach my $mirna (@{$chromosomes{$chr}{$strand}})
	    {
		my $start = ($mirna->{start}-$clustersize > 0) ? $mirna->{start}-$clustersize : 0;
		my $stop  = $mirna->{stop}+$clustersize;
		for(my $i=$start; $i<=$stop; $i++)
		{
		    $dna_strand[$i] = "c";
		}
	    }

	    # find cluster borders
	    my ($cluster_start, $cluster_stop);
	    my @cluster = ();

	    for (my $i=0; $i<@dna_strand; $i++)
	    {
		next unless (defined $dna_strand[$i]);
		if ($dna_strand[$i] eq "c")
		{
		    $cluster_start=$i;
		    for (my $j=$cluster_start+1; $j<@dna_strand; $j++)
		    {
			if (! defined $dna_strand[$j] || $dna_strand[$j] ne "c")
			{
			    $cluster_stop=$j-1;
			    $i = $j;
			    last;
			}

		    }

		    unless (defined $cluster_stop)
		    {
			$cluster_stop = @dna_strand-1;
			$i = $cluster_stop;
		    }

		    push(@cluster, {start => $cluster_start, stop => $cluster_stop});
		    ($cluster_start, $cluster_stop) = (undef, undef);
		}
	    }

	    # cluster are extracted... Find all miRNAs inside the cluster
	    foreach my $clust (@cluster)
	    {
		my @clustered_miRNAs = grep {$_->{start}>=$clust->{start} && $_->{stop} <= $clust->{stop}} (@{$chromosomes{$chr}{$strand}});
		if (@clustered_miRNAs == 1)
		{
		    print STDERR "Single miRNA cluster for miRNA $clustered_miRNAs[0]{mirna} on chromosome $chr($strand)\n";
		} else {
		    printf "Cluster of miRNAs (border-size: %d) on chromosome %s(%s) (%d-%d) containing %d miRNAs: %s\n", $clustersize, $chr, $strand, $clust->{start}, $clust->{stop}, (@clustered_miRNAs+0), join(", ", map {$_->{mirna}} (@clustered_miRNAs));
		}
	    }

	} else {
	    die "Should never happen";
	}
    }
}
