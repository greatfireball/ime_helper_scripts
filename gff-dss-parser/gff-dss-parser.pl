#!/usr/bin/env perl

use strict;
use warnings;

use Bio::Tools::GFF;

my $file = "manduca_genomic_regions.gff";

my $gffio = Bio::Tools::GFF->new(-file => $file);

