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

GetOptions(
    "mascot=s" => \$mascot_id_file
    ) || die;

my %mascot_ids = ();

open(FH, "<", $mascot_id_file) || die "Unable to open Mascot ID File '$mascot_id_file': $!\n";
while(<FH>)
{
    chomp;
    if ($_ =~ /^.*?(\d+)\.(p\d+)/)
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
    
