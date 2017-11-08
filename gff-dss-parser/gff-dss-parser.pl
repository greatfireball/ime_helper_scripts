#!/usr/bin/env perl

use strict;
use warnings;

use Text::CSV;
use Bio::Range;

my %gff = ();
my $counter = 0;

my $file;
my @csv;
my @new;

use Getopt::Long;
use Storable qw(dclone);

GetOptions( 'gff=s' => \$file,
	    'csv=s@' => \@csv,
	    'out=s@' => \@new
    ) || die;

@csv = split(",", join(",", @csv));
@new = split(",", join(",", @new));

my $toolname = "reannotator";

unless ($file && -e $file)
{
    die "Missing input gff. Use --gff parameter!\n";
}

unless (@csv)
{
    die "Missing input cvs file(s). Use --cvs parameter!\n";
}

unless (@new)
{
    warn "Missing output file(s). Deriving output names from input csv files. Otherwise use --out parameter!\n";
    @new = map {$_.".out"} (@csv);
}

unless (@csv == @new)
{
    die "Number of input and output files differ!\n";
}

foreach my $csv (@csv)
{
    unless ($csv && -e $csv)
    {
	die "Missing input cvs file '$csv'!\n";
    }
}

foreach my $new (@new)
{
    unless (! -e $new)
    {
	die "Output file '$new' exists. Please use other name using --out parameter!\n";
    }
}

open(FH, "<", $file) || die "$!\n";

my %sort_types = (
    CDS => 50,
    downstream => 5,
    exon => 40,
    five_prime_UTR => 45,
    gene => 10,
    'gene-intron' => 15,
    intergenic => 2,
    intron => 25,
    mRNA => 20,
    three_prime_UTR => 45,
    upstream => 5,
    non_canonical_five_prime_splice_site => 1,
    non_canonical_three_prime_splice_site => 1,
    transcript => 1
    );

sub sort_gff_types
{
    my $a_cleaned = $a;
    $a_cleaned =~ s/[.+-]$//;
    my $b_cleaned = $b;
    $b_cleaned =~ s/[.+-]$//;

    $sort_types{$b_cleaned} <=> $sort_types{$a_cleaned} || $a cmp $b;
}

warn "Sorting in the following order (top-->down): ", join(", ", sort sort_gff_types (keys %sort_types)), "\n";

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

warn "Generating UTR annotations\n";

generate_UTR_annotation(\%gff);

warn "Finished\n";

warn "Generating intron annotations\n";

generate_intron_annotation(\%gff);

warn "Finished\n";

warn "Sorting information\n";

sort_gff(\%gff);

warn "Finished\n";

print_gff(\%gff);

for(my $i=0; $i<@csv; $i++)
{
    my $csv = $csv[$i];
    my $new = $new[$i];

    warn sprintf("##### Starting work on input file '%s' for output file '%s'\n", $csv, $new);

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
		$types{$_->{orig}{type}.$_->{orig}{strand}}++;
	    }

	    my @sorted_types = sort sort_gff_types (keys %types);

	    # check if no annoation was returned
	    unless (@sorted_types)
	    {
		@sorted_types = ("none");
		$types{"none"} = 1;
	    }

	    push(@{$row}, $sorted_types[0], join(":", map { sprintf("%s(%d)", $_, $types{$_}) } (@sorted_types)));
	} else {
	    # no entry was found for chromosome

	    warn "Missing entry for chromosome '$chr'... Assume no annoation\n";

	    push(@{$row}, "chr_not_annotated", "chr_not_annotated");
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
}

sub get_parent_child_relationship
{
    my $gff = shift;

    my %relationship = ();

    foreach my $chr (keys %{$gff})
    {
	for(my $i=0; $i < @{$gff->{$chr}}; $i++)
	{
	    my $entry = $gff->{$chr}[$i];

	    if (exists $entry->{attributes}{ID} && @{$entry->{attributes}{ID}}==1)
	    {
		my $parent = $entry->{attributes}{ID}[0];
		$relationship{$chr}{$parent}{parent} = $i;
	    }

	    if (exists $entry->{attributes}{Parent})
	    {
		foreach my $parent (@{$entry->{attributes}{Parent}})
		{
		    push(@{$relationship{$chr}{$parent}{children}}, $i);
		}
	    }
	}

	# # deleting entries without children
	# my $chromosome_set = $relationship{$chr};
	# foreach my $entry2delete (grep {! exists $chromosome_set->{$_}{children}} (keys %{$chromosome_set}))
	# {
	#     delete $chromosome_set->{$entry2delete};
	# }
    }

    return \%relationship;
}

sub generate_UTR_annotation
{
    my $gff = shift;

    my $relationship = get_parent_child_relationship($gff);;

    foreach my $chr (keys %{$relationship})
    {
	foreach my $parent (keys %{$relationship->{$chr}})
	{
	    # check if the group contains CDS and exon information
	    my @CDS = sort {
		$a->{start} <=> $b->{start} ||
		    $a->{stop} <=> $b->{stop}
	    } grep {
		$_->{type} eq "CDS"
	    } map {{%{$gff->{$chr}[$_]}}} (@{$relationship->{$chr}{$parent}{children}});

	    my @exons = sort {
		$a->{start} <=> $b->{start} ||
		    $a->{stop} <=> $b->{stop}
	    } grep { $_->{type} eq "exon" } map {{%{$gff->{$chr}[$_]}}} (@{$relationship->{$chr}{$parent}{children}});

	    next unless (@exons && @CDS);

	    my @new_entries = ();
	    # search for first exon containing CDS
	    my $first_CDS = $CDS[0];
	    my $reverse = 0;
	    if ($first_CDS->{strand} eq "-")
	    {
		$reverse = 1;
	    }
	    foreach my $exon (@exons)
	    {
		my $new_entry = dclone($exon);
		$new_entry->{type} = ($reverse) ? "three_prime_UTR" : "five_prime_UTR";
		$new_entry->{source} = $toolname;
		delete $new_entry->{attributes}{Name};
		delete $new_entry->{attributes}{ID};

		if ($first_CDS->{start} >= $exon->{start} && $first_CDS->{start} <= $exon->{stop})
		{
		    if ($first_CDS->{start} > $exon->{start})
		    {
			# found a match, therefore partial UTR
			$new_entry->{stop} = $first_CDS->{start}-1;
			push(@new_entries, $new_entry);
		    }

		    last;
		} else {
		    # no match, therefore complete UTR exon
		    push(@new_entries, $new_entry);
		}
	    }

	    my $last_CDS = $CDS[@CDS-1];
	    foreach my $exon (reverse @exons)
	    {
		my $new_entry = dclone($exon);
		$new_entry->{type} = ($reverse) ? "five_prime_UTR" : "three_prime_UTR";
		$new_entry->{source} = $toolname;
		delete $new_entry->{attributes}{Name};
		delete $new_entry->{attributes}{ID};

		if ($last_CDS->{stop} >= $exon->{start} && $last_CDS->{stop} <= $exon->{stop})
		{
		    if ($last_CDS->{stop} < $exon->{stop})
		    {
			# found a match, therefore partial UTR
			$new_entry->{start} = $last_CDS->{stop}+1;
			push(@new_entries, $new_entry);
		    }

		    last;
		} else {
		    # no match, therefore complete UTR exon
		    push(@new_entries, $new_entry);
		}
	    }

	    push(@{$gff->{$chr}}, @new_entries);
	}
    }
    return 1;
}

sub generate_intron_annotation
{
    my $gff = shift;

    my $relationship = get_parent_child_relationship($gff);;

    foreach my $chr (keys %{$relationship})
    {
	foreach my $parent (keys %{$relationship->{$chr}})
	{
	    # check if the group contains exon information
	    my @exons = sort {
		$a->{start} <=> $b->{start} ||
		    $a->{stop} <=> $b->{stop}
	    } grep { $_->{type} eq "exon" } map {{%{$gff->{$chr}[$_]}}} (@{$relationship->{$chr}{$parent}{children}});

	    next unless (@exons>1);

	    my @introns = ();

	    # search for first exon containing CDS
	    for(my $i=0; $i<@exons-1; $i++)
	    {
		my $intron = dclone($exons[$i]);
		$intron->{type} = "intron";
		$intron->{source} = $toolname;
		delete $intron->{attributes}{Name};
		delete $intron->{attributes}{ID};

		$intron->{start}++;
		$intron->{stop} = $exons[$i+1]{start}-1;

		push(@introns, $intron);
	    }

	    push(@{$gff->{$chr}}, @introns);
	}
    }
    return 1;
}

sub sort_gff
{
    my $gff = shift;

    my $parent_child_rel = get_parent_child_relationship($gff);

    foreach my $chr (keys %{$gff})
    {
	@{$gff->{$chr}} = map {
	    {
		range => Bio::Range->new(-start => $_->{start}, -end => $_->{stop}, -strand => 0),
		orig => $_
	    }
	} sort {
	    $sort_types{$a->{type}} <=> $sort_types{$b->{type}} ||
	$a->{start} <=> $b->{start} ||
	$a->{stop} <=> $b->{stop}
	} @{$gff->{$chr}};
    }
}

sub print_gff
{
    my $gff = shift;

    foreach my $chr (sort keys %{$gff})
    {
	foreach my $entry (@{$gff->{$chr}})
	{
	    print join("\t", (map {$entry->{orig}{$_}} (qw(chromosome source type start stop score strand phase))), join(";", map { sprintf("%s=%s", $_, join(",", @{$entry->{orig}{attributes}{$_}})) } (keys %{$entry->{orig}{attributes}}))), "\n";
	}
    }
}
