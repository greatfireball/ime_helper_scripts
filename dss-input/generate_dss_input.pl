#!/usr/bin/env perl
use strict;
use warnings;

use IO::Uncompress::Gunzip qw($GunzipError);

=head1 fastq-shuffle.pl

A small program to convert and combine bismark coverage output files
into DSS input files.

=head1 SYNOPSIS

    generate_dss_input bismark1.cov.gz bismark2.cov

=head1 OUTPUT

The output is a csv file containing a header and the columns

=over 4

=item 1. chr

=item 2. pos

=item 3. N

=item 4. X

=back

as required by DSS. The counts for total coverage and methylated
coverage of all input files is combined and the sum will be printed to
stdout.

=head1 SEE ALSO

=over 4

=item Bismark

A mapper for bs-seq C<>

=item DSS

A R package for differential methylation experiments C<https://www.bioconductor.org/packages/release/bioc/html/DSS.html>

=back

=head1 AUTHOR

This script was written by Frank FE<ouml>rster C<frank.foerster@ime.fraunhofer.de>.

=head1 CHANGELOG

=over 4

=item 2017-10-16 v0.1

First working version.

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

use version 0.77; our $VERSION = version->declare("v0.1.0");

my @infiles = @ARGV;
my %dat     = ();

foreach my $file (@infiles) {
    my $z;
    if ($file =~ /\.gz$/)
    {
	$z = new IO::Uncompress::Gunzip $file || die "gunzip files $GunzipError\n";
    } else {
	open( $z, "<", $file ) || die "$!";
    }
    while (<$z>) {
        chomp $_;
        my ( $chr, $start, $end, $methyl_percent, $count_methyl,
            $count_notmethyl )
          = split( /\t/, $_ );

        for ( my $pos = $start ; $pos <= $end ; $pos++ ) {
            if ( !exists $dat{$chr}{$pos}{$file} )
            {
                $dat{$chr}{$pos}{$file} =
                  [ $count_methyl + $count_notmethyl, $count_methyl ];
            }
            else {
                die "Something wrong for $chr $pos $file";
            }
        }
    }
    close($z) || die "$!";
}

# printed overall sums sorted by chromosome followed by position
print join( "\t", ( "chr", "pos", "N", "X" ) ), "\n";
foreach my $chr ( sort keys %dat ) {
    foreach my $pos ( sort { $a <=> $b } ( keys %{ $dat{$chr} } ) ) {
        my ( $total_count, $methyl_count ) = ( 0, 0 );
        foreach my $file (@infiles) {
            if ( exists $dat{$chr}{$pos}{$file} ) {
                $total_count  += $dat{$chr}{$pos}{$file}[0];
                $methyl_count += $dat{$chr}{$pos}{$file}[1];
            }
        }

        print join( "\t", ( $chr, $pos, $total_count, $methyl_count ) ), "\n";
    }
}
