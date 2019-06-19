#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

use Log::Log4perl;

# Configuration in a string ...
my $conf = q(
  log4perl.category.Foo.Bar          = INFO, Logfile, Screen

  log4perl.appender.Logfile          = Log::Log4perl::Appender::File
  log4perl.appender.Logfile.filename = test.log
  log4perl.appender.Logfile.layout   = Log::Log4perl::Layout::PatternLayout
  log4perl.appender.Logfile.layout.ConversionPattern = [%r] %F %L %m%n

  log4perl.appender.Screen         = Log::Log4perl::Appender::Screen
  log4perl.appender.Screen.stderr  = 0
  log4perl.appender.Screen.layout = Log::Log4perl::Layout::SimpleLayout
);

# ... passed as a reference to init()
Log::Log4perl::init( \$conf );
my $log = Log::Log4perl::get_logger("Foo::Bar");

my $mascot_id_file = "";
my $interproscan_file = "";
my $toxprot_file      = "";

GetOptions(
    "mascot=s"      => \$mascot_id_file,
    "interpro=s"    => \$interproscan_file,
    "toxprot=s"     => \$toxprot_file
    ) || die;

my %mascot_ids = ();

open(FH, "<", $mascot_id_file) || die "Unable to open Mascot ID File '$mascot_id_file': $!\n";
while(<FH>)
{
    chomp;
    if ($_ =~ /^\S*?(\d+)\.(p\d+)/)
    {
	my ($id, $protein) = ($1, $2);
	push(@{$mascot_ids{$id}{orig}}, $_);
	$mascot_ids{$id}{proteins}{$protein}++;
    } else {
	die "Should never happen for Mascot ID File line $. : '$_'\n";
    }
}
close(FH) || die "Unable to close Mascot ID File '$mascot_id_file': $!\n";

$log->info(sprintf("Found %d ids from mascot run", (keys %mascot_ids)+0));

open(FH, "<", $interproscan_file) || die "Unable to open Interproscan file '$interproscan_file': $!\n";
while(<FH>)
{
    if ($_ =~ /^\S*?(\d+)\.(p\d+)/)
    {
	my ($id, $protein) = ($1, $2);
	if (exists $mascot_ids{$id})
	{
	    push(@{$mascot_ids{$id}{interpro}}, $_)
	}
    }
}
close(FH) || die "Unable to close Interproscan file '$interproscan_file': $!\n";

my @ids_with_interpro = ();
foreach my $id (keys %mascot_ids)
{
    if (exists $mascot_ids{$id}{interpro})
    {
	push(@ids_with_interpro, $id);
    } else {
	$mascot_ids{$id}{interpro} = [];
    }
}
my $num_interpros = @ids_with_interpro+0;
my $percent_interpro = $num_interpros/((keys %mascot_ids)+0)*100;
$log->info(sprintf("Found %d entries with interproscan information (%.1f %%)", $num_interpros, $percent_interpro));

open(FH, "<", $toxprot_file) || die "Unable to open toxprot file '$toxprot_file': $!\n";
my @fieldnames = ();
while(<FH>)
{
    if (/^#/)
    {
	chomp;
	$_ =~ s/^#//;
	@fieldnames = split(/\t/, $_);
	next;
    }
    chomp;
    # example line
    # #query  query_lengh     subject subject_length  identity        similarity      query_start     query_end       subject_start   subject_end     query_coverage  subject_coverage        bitscore        description
    # pitu_v1_bf_011616       1837    sp|J3SDX8|LICH_CROAD    400     36.412  57.18   332     1459    26      397     61.35   92.75   225     Putative lysosomal acid lipase/cholesteryl ester hydrolase OS=Crotalus adamanteus OX=8729 PE=2 SV=1
    my @fields = split(/\t/, $_);
    my %data = ();
    @data{@fieldnames} = @fields;

    if ($data{query} =~ /^\S+_(\d+)$/)
    {
	my ($id) = ($1);
	if (exists $mascot_ids{$id})
	{
	    # update if it is a better hit or no toxprot hit was stored
	    if ( (! exists $mascot_ids{$id}{toxprot}) || ($mascot_ids{$id}{toxprot}{bitscore} < $data{bitscore}))
	    {
		$mascot_ids{$id}{toxprot} = \%data;
	    }
	}
    }
}
close(FH) || die "Unable to close toxprot file '$toxprot_file': $!\n";

my @ids_with_toxprot = ();
foreach my $id (keys %mascot_ids)
{
    if (exists $mascot_ids{$id}{toxprot})
    {
	push(@ids_with_toxprot, $id);
    } else {
	$mascot_ids{$id}{toxprot} = {};
    }
}
my $num_toxprot = @ids_with_toxprot+0;
my $percent_toxprot = $num_toxprot/((keys %mascot_ids)+0)*100;
$log->info(sprintf("Found %d entries with toxprot information (%.1f %%)", $num_toxprot, $percent_toxprot));
