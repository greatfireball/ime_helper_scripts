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

GetOptions(
    "mascot=s"      => \$mascot_id_file,
    "interpro=s"    => \$interproscan_file
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
