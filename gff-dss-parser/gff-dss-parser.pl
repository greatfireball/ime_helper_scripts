#!/usr/bin/env perl

use strict;
use warnings;

use Text::CSV;
use Bio::Range;

my %gff = ();
my $counter = 0;

my $file = "manduca_genomic_regions.gff";
my $csv = "dmrs_sm0.05.csv";
my $new = "dmrs_sm0.05.new.csv";

open(FH, "<", $file) || die "$!\n";

my %sort_types = (
    CDS => 50,
    downstream => 5,
    exon => 40,
    five_prime_UTR => 30,
    gene => 10,
    'gene-intron' => 15,
    intergenic => 5,
    intron => 25,
    mRNA => 20,
    three_prime_UTR => 30,
    upstream => 5
    );

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

    push(@{$gff{$chr}}, {
	chromosome => $chr,
	source => $source,
	type => $type,
	start => $start,
	stop => $stop,
	score => $score,
	strand => $strand,
	phase => $phase,
	attributes => \%attributes
	 });
}

close(FH) || die "$!\n";

warn "Sorting information\n";

foreach my $chr (keys %gff)
{
    @{$gff{$chr}} = map {{range => Bio::Range->new(-start => $_->{start}, -end => $_->{stop}, -strand => 0), orig => $_ }} sort {$a->{start} <=> $b->{start} || $a->{stop} <=> $b->{stop} || $a->{type} cmp $b->{type}} @{$gff{$chr}};
}

warn "Finished\n";

my @rows = ();
my $csv_parser = Text::CSV->new ( { binary => 1, sep_char => "\t" } )  # should set binary attribute.
    or die "Cannot use CSV: ".Text::CSV->error_diag ();

open(my $fh, "<:encoding(utf8)", $csv) or die "'$csv': $!";

my $col_names = $csv_parser->getline( $fh );
$csv_parser->column_names($col_names);

while ( my $row = $csv_parser->getline( $fh ) )
{
    # correct sequence name
    if ($row->[0] =~ /^gi\|\d+\|gb\|([^|]+)\|$/)
    {
	$row->[0] = $1;
    }

    # are there annotations for the given region
    my $chr = $row->[0];
    my $start = $row->[1]+0;
    my $stop = $row->[2]+0;

    my $target = Bio::Range->new(-start => $start, -end => $stop, -strand => 0);

    if (exists $gff{$chr})
    {
	# find annotations
	my @annotations = grep { my $r = $_->{range}; $r->overlaps($target); } (@{$gff{$chr}});

	my %types = ();
	foreach (@annotations)
	{
	    $types{$_->{orig}{type}}++;
	}

	my @sorted_types = sort { $sort_types{$b} <=> $sort_types{$a} || $a cmp $b } (keys %types);

	push(@{$row}, $sorted_types[0], join(":", map { sprintf("%s(%d)", $_, $types{$_}) } (@sorted_types)));
    } else {
	# no entry was found for chromosome

	warn "Missing entry for chromosome '$chr'... Assume no annoation\n";

	push(@{$row}, "none", "none");
    }

    push(@rows, $row);
}
$csv_parser->eof or $csv_parser->error_diag();
close($fh) || die "'$csv': $!";

# write output
# add a field to the col_names
$col_names = [ @{$col_names}, "annotation", "all_annotations" ];
$csv_parser->eol ("\n");
open($fh, ">:encoding(utf8)", $new) or die "'$new': $!";
$csv_parser->print ($fh, $col_names);
$csv_parser->print ($fh, $_) for @rows;
close($fh) or die "'$new': $!";
