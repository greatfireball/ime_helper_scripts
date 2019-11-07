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
	die "Epic fail\n" if (exists $id2seq->{$seq->id} && $id2seq->{$seq->id} ne $seq->seq);

	$id2seq->{$seq->id} = $seq->seq;
	push(@{$mapping->{$seq->seq}}, $seq->id);
    }
}

open(FH, "<", $mascot_id_file) || die "Unable to open file '$mascot_id_file': $!\n";
my $i=0;
while(<FH>)
{
    chomp;
    warn "Missing sequence for $_\n" unless (exists $id2seq->{$_});

    $i++;
    foreach my $id (sort @{$mapping->{$id2seq->{$_}}})
    {
	printf "Group_%06d\t%s\n", $i, $id;
    }
}
close(FH) || die "Unable to close file '$mascot_id_file: $!\n";

# my @seqs = sort { length($b) <=> length($a) || $a cmp $b } (keys (%{$mapping}));
# for(my $i=0; $i<@seqs; $i++)
# {
#     foreach my $id (sort @{$mapping->{$seqs[$i]}})
#     {
# 	printf "Group_%06d\t%s\n", $i, $id;
#     }
# }
