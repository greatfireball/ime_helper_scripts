#!/usr/bin/env perl
use strict;
use warnings;

use IO::Uncompress::Gunzip qw($GunzipError);

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
