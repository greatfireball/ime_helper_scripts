#!/usr/bin/env perl

use strict;
use warnings;

use Bio::SeqIO;

use Getopt::Long;

my @seqfiles       = ();
my $mascot_id_file = "";

GetOptions("sequences=s@" => \@seqfiles,
	   "mascot=s"     => \$mascot_id_file) || die "Error in command line arguments\n";

@seqfiles = split(",", join(",", @seqfiles));

my $mapping = {};
my $id2seq  = {};

foreach my $seqfile (@seqfiles)
{
    my $inseq = Bio::SeqIO->new(-file => $seqfile);
    while (my $seq = $inseq->next_seq)
    {
	die "Epic fail\n" if (exists $id2seq->{$seq->id} && $id2seq->{$seq->is} ne $seq->seq);

	$id2seq->{$seq->id} = $seq->seq;
	push(@{$mapping->{$seq->seq}}, $seq->id);
    }
}
