#!/usr/bin/env perl

use strict;
use warnings;

use Digest::MD5;

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

use version 0.77; our $VERSION = version->declare('v1.2.0');

# ... passed as a reference to init()
Log::Log4perl::init( \$conf );
my $log = Log::Log4perl::get_logger("Foo::Bar");

my $mascot_id_file      = "";
my $interproscan_file   = "";
my $toxprot_file        = "";
my $ncbi_file           = "";
my $quantification_file = "";
my $hmmer_file          = "";
my $jackhmmer_file      = "";
my $deseq_normalized    = "";
my @deseq_results       = ();
my @manual_annotations  = ();
my $interpro_out        = "interpro.tsv";
my $table_out           = "output.tsv";

GetOptions(
    "mascot=s"         => \$mascot_id_file,
    "interpro=s"       => \$interproscan_file,
    "toxprot=s"        => \$toxprot_file,
    "ncbi=s"           => \$ncbi_file,
    "quantification=s" => \$quantification_file,
    "hmmer=s"          => \$hmmer_file,
    "jackhmmer=s"      => \$jackhmmer_file,
    "deseqnormal=s"    => \$deseq_normalized,
    "deseqdiffseq=s@"  => \@deseq_results,
    "manual=s@"        => \@manual_annotations,
    "out=s"            => \$table_out,
    "out_interpro=s"   => \$interpro_out
    ) || die;

@manual_annotations = split(/,/, join(",", @manual_annotations));
@deseq_results = split(/,/, join(",", @deseq_results));

my %mascot_ids = ();

my $ctx = Digest::MD5->new;
open(FH, "<", $mascot_id_file) || die "Unable to open Mascot ID File '$mascot_id_file': $!\n";
while(<FH>)
{
    $ctx->add($_);
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

$log->info(sprintf("Found %d ids from mascot run (MD5-Input: %s)", (keys %mascot_ids)+0, $ctx->hexdigest));

for my $manual (@manual_annotations)
{
    $ctx = Digest::MD5->new;
    open(FH, "<", $manual) || die "Unable to open manual annotation file '$manual': $!\n";
    while(<FH>)
    {
	$ctx->add($_);
	if ($_ =~ /^\S+_(\d+)\t+(.+)/)
	{
	    my ($id, $annotation) = ($1, $2);
	    if (exists $mascot_ids{$id})
	    {
		push(@{$mascot_ids{$id}{manual}}, $annotation)
	    }
	}
    }
    close(FH) || die "Unable to close manual annotation file '$manual': $!\n";

    my @ids_with_manual = ();
    foreach my $id (keys %mascot_ids)
    {
	if (exists $mascot_ids{$id}{manual})
	{
	    push(@ids_with_manual, $id);
	} else {
	    $mascot_ids{$id}{manual} = [];
	}
    }
    my $num_manuals = @ids_with_manual+0;
    my $percent_manual = $num_manuals/((keys %mascot_ids)+0)*100;
    $log->info(sprintf("Found %d entries with manual annoation from file '%s' (%.1f %%) (MD5-Input: %s)", $num_manuals, $manual, $percent_manual, $ctx->hexdigest));
}

$ctx = Digest::MD5->new;
open(FH, "<", $interproscan_file) || die "Unable to open Interproscan file '$interproscan_file': $!\n";
while(<FH>)
{
    $ctx->add($_);
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
$log->info(sprintf("Found %d entries with interproscan information (%.1f %%) (MD5-Input: %s)", $num_interpros, $percent_interpro, $ctx->hexdigest));

$ctx = Digest::MD5->new;
open(FH, "<", $hmmer_file) || die "Unable to open hmmer file '$hmmer_file': $!\n";
while(<FH>)
{
    $ctx->add($_);

    next if (/^#/);

    chomp;
    my @fields = split(/\s+/, $_);

    if ($fields[0] =~ /^\S*?(\d+)\.(p\d+)/)
    {
	my ($id) = ($1);
	if (exists $mascot_ids{$id})
	{
	    my $bitscore = $fields[5];
	    my $name     = $fields[2];

	    next if ($bitscore < 15);

	    if ( (! exists $mascot_ids{$id}{hmmer}) || $mascot_ids{$id}{hmmer}{bitscore} < $bitscore)
	    {
		$mascot_ids{$id}{hmmer}= { bitscore => $bitscore, hmmerhit => $name };
	    }
	}
    }
}
close(FH) || die "Unable to close hmmer file '$hmmer_file': $!\n";

my @ids_with_hmmer = ();
foreach my $id (keys %mascot_ids)
{
    if (exists $mascot_ids{$id}{hmmer})
    {
	push(@ids_with_hmmer, $id);
    } else {
	$mascot_ids{$id}{hmmer} = {};
    }
}
my $num_hmmer = @ids_with_hmmer+0;
my $percent_hmmer = $num_hmmer/((keys %mascot_ids)+0)*100;
$log->info(sprintf("Found %d entries with hmmer information (%.1f %%) (MD5-Input: %s)", $num_hmmer, $percent_hmmer, $ctx->hexdigest));

$ctx = Digest::MD5->new;
open(FH, "<", $jackhmmer_file) || die "Unable to open jackhmmer file '$jackhmmer_file': $!\n";
while(<FH>)
{
    $ctx->add($_);

    next if (/^#/);

    chomp;
    my @fields = split(/\s+/, $_);

    if ($fields[0] =~ /^\S*?(\d+)\.(p\d+)/)
    {
	my ($id) = ($1);
	if (exists $mascot_ids{$id})
	{
	    my $bitscore = $fields[5];
	    my $name     = $fields[2];

	    next if ($bitscore < 15);

	    if ( (! exists $mascot_ids{$id}{jackhmmer}) || $mascot_ids{$id}{jackhmmer}{bitscore} < $bitscore)
	    {
		$mascot_ids{$id}{jackhmmer}= { bitscore => $bitscore, jackhmmerhit => $name };
	    }
	}
    }
}
close(FH) || die "Unable to close jackhmmer file '$jackhmmer_file': $!\n";

my @ids_with_jackhmmer = ();
foreach my $id (keys %mascot_ids)
{
    if (exists $mascot_ids{$id}{jackhmmer})
    {
	push(@ids_with_jackhmmer, $id);
    } else {
	$mascot_ids{$id}{jackhmmer} = {};
    }
}
my $num_jackhmmer = @ids_with_jackhmmer+0;
my $percent_jackhmmer = $num_jackhmmer/((keys %mascot_ids)+0)*100;
$log->info(sprintf("Found %d entries with jackhmmer information (%.1f %%) (MD5-Input: %s)", $num_jackhmmer, $percent_jackhmmer, $ctx->hexdigest));

$ctx = Digest::MD5->new;
open(FH, "<", $quantification_file) || die "Unable to open quantification file '$quantification_file': $!\n";
my @fieldnames = ();
while(<FH>)
{
    $ctx->add($_);

    if (/^#/)
    {
	chomp;
	$_ =~ s/^#//;
	@fieldnames = split(/\t/, $_);
	next;
    }
    chomp;
    my @fields = split(/\t/, $_);
    my %data = ();
    @data{@fieldnames} = @fields;

    if ($data{transcript} =~ /^\S+_(\d+)$/)
    {
	my ($id) = ($1);
	if (exists $mascot_ids{$id})
	{
	    die "Double quantification result for $id\n" if (exists $mascot_ids{$id}{quantification});
	    $mascot_ids{$id}{quantification} = \%data;
	}
    }
}
close(FH) || die "Unable to close quantification file '$quantification_file': $!\n";

my @ids_with_quantification = ();
foreach my $id (keys %mascot_ids)
{
    if (exists $mascot_ids{$id}{quantification})
    {
	push(@ids_with_quantification, $id);
    } else {
	$mascot_ids{$id}{quantification} = {};
    }
}
my $num_quantification = @ids_with_quantification+0;
my $percent_quantification = $num_quantification/((keys %mascot_ids)+0)*100;
$log->info(sprintf("Found %d entries with quantification information (%.1f %%) (MD5-Input: %s)", $num_quantification, $percent_quantification, $ctx->hexdigest));

$ctx = Digest::MD5->new;
open(FH, "<", $toxprot_file) || die "Unable to open toxprot file '$toxprot_file': $!\n";
while(<FH>)
{
    $ctx->add($_);

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
$log->info(sprintf("Found %d entries with toxprot information (%.1f %%) (MD5-Input: %s)", $num_toxprot, $percent_toxprot, $ctx->hexdigest));

$ctx = Digest::MD5->new;
open(FH, "<", $ncbi_file) || die "Unable to open ncbi file '$ncbi_file': $!\n";
while(<FH>)
{
    $ctx->add($_);
    if (/^#/)
    {
	chomp;
	$_ =~ s/^#//;
	@fieldnames = split(/\t/, $_);
	next;
    }
    chomp;
    my @fields = split(/\t/, $_);
    my %data = ();
    @data{@fieldnames} = @fields;

    if ($data{query} =~ /^\S+_(\d+)$/)
    {
	my ($id) = ($1);
	if (exists $mascot_ids{$id})
	{
	    # update if it is a better hit or no toxprot hit was stored
	    if ( (! exists $mascot_ids{$id}{ncbi}) || ($mascot_ids{$id}{ncbi}{bitscore} < $data{bitscore}))
	    {
		$mascot_ids{$id}{ncbi} = \%data;
	    }
	}
    }
}
close(FH) || die "Unable to close ncbi file '$ncbi_file': $!\n";

my @ids_with_ncbi = ();
foreach my $id (keys %mascot_ids)
{
    if (exists $mascot_ids{$id}{ncbi})
    {
	push(@ids_with_ncbi, $id);
    } else {
	$mascot_ids{$id}{ncbi} = {};
    }
}
my $num_ncbi = @ids_with_ncbi+0;
my $percent_ncbi = $num_ncbi/((keys %mascot_ids)+0)*100;
$log->info(sprintf("Found %d entries with ncbi information (%.1f %%) (MD5-Input: %s)", $num_ncbi, $percent_ncbi, $ctx->hexdigest));

$ctx = Digest::MD5->new;
@fieldnames = ();
open(FH, "<", $deseq_normalized) || die "Unable to open DESeq normalized file '$deseq_normalized': $!\n";
while(<FH>)
{
    $ctx->add($_);

    if ($. == 1)
    {
	chomp;
	@fieldnames = split(/\t/, $_);
	# add an ID field
	unshift(@fieldnames, "id");
	next;
    }
    chomp;
    my @fields = split(/\t/, $_);
    my %data = ();
    @data{@fieldnames} = @fields;

    if ($data{id} =~ /^\S+_(\d+)$/)
    {
	my ($id) = ($1);
	if (exists $mascot_ids{$id})
	{
	    $mascot_ids{$id}{deseqnormalized} = \%data;
	}
    }
}
close(FH) || die "Unable to close DESeq normalized file '$deseq_normalized': $!\n";

my @ids_with_normalized_counts = ();
foreach my $id (keys %mascot_ids)
{
    if (exists $mascot_ids{$id}{deseqnormalized})
    {
	push(@ids_with_normalized_counts, $id);
    } else {
	$mascot_ids{$id}{deseqnormalized} = {
	    id         => $id,
	    venomgland => "-",
	    femalebody => "-",
	    malebody   => "-"
	};
    }
}
my $num_deseqnormalized = @ids_with_normalized_counts+0;
my $percent_deseqnormalized = $num_deseqnormalized/((keys %mascot_ids)+0)*100;
$log->info(sprintf("Found %d entries with DESeq normalized information (%.1f %%) (MD5-Input: %s)", $num_deseqnormalized, $percent_deseqnormalized, $ctx->hexdigest));

foreach my $diffseq_file (@deseq_results)
{
    @fieldnames = ();
    my ($experiment) = $diffseq_file =~ /([^_]+)\.mat$/;

    $ctx = Digest::MD5->new;
    open(FH, "<", $diffseq_file) || die "Unable to open DESeq diffseq file '$diffseq_file': $!\n";
    while(<FH>)
    {
	$ctx->add($_);
	if ($. == 1)
	{
	    chomp;
	    @fieldnames = split(",", $_);
	    # skip first field
	    shift @fieldnames;
	    next;
	}
	chomp;
	my @fields = split(",", $_);
	shift @fields;
	my %data = ();
	@data{@fieldnames} = @fields;

	if ($data{id} =~ /^\S+_(\d+)$/)
	{
	    my ($id) = ($1);
	    if (exists $mascot_ids{$id})
	    {
		$mascot_ids{$id}{diffseq}{$experiment} = {
		    id => $id,
		    'log2foldchange' => $data{log2FoldChange},
		    'significant'    => (($data{padj} <= 0.05) ? "YES" : "NO"),
		    'padj'           => $data{padj}
		};
	    }
	}
    }
    close(FH) || die "Unable to close DESeq diffseq '$diffseq_file': $!\n";

    my @ids_with_diffseq_counts = ();
    foreach my $id (keys %mascot_ids)
    {
	if (exists $mascot_ids{$id}{diffseq}{$experiment})
	{
	    push(@ids_with_diffseq_counts, $id);
	} else {
	    $mascot_ids{$id}{diffseq}{$experiment} = {
		id               => $id,
		'log2foldchange' => "-",
		'significant'    => "NO",
		'padj'           => "-"
	    };
	}
    }
    my $num_diffseq = @ids_with_diffseq_counts+0;
    my $percent_diffseq = $num_diffseq/((keys %mascot_ids)+0)*100;
    $log->info(sprintf("Found %d entries with DESeq diffseq information (%.1f %%) for input file '%s' (MD5-Input: %s)", $num_diffseq, $percent_diffseq, $diffseq_file, $ctx->hexdigest));
}

my $ctx_interpro = Digest::MD5->new;
my $ctx_output   = Digest::MD5->new;
open(FH, ">", $table_out) || die "Unable to open output table file '$table_out': $!\n";
open(INTERPRO, ">", $interpro_out) || die "Unable to open interpro output table file '$interpro_out': $!\n";

my @diffseq_fields = ();
foreach my $diffseq_file (@deseq_results)
{
    my ($experiment) = $diffseq_file =~ /([^_]+)\.mat$/;

    push(@diffseq_fields, map { $experiment."_".$_ } (qw(log2foldchange significant padj)));
}

my $line = join("\t",
	      "transcript",
	      "vg_cov",
	      "vg_tpm",
	      "bf_cov",
	      "bf_tpm",
	      "bm_cov",
	      "bm_tpm",
	      "vg_DESeq_normalized",
	      "bf_DESeq_normalized",
	      "bm_DESeq_normalized",
	      @diffseq_fields,
	      "Venom Related Protein Description",
	      "ToxProt_Subject",
	      "ToxProt_Description",
	      "ToxProt_Bitscore",
	      "NCBI_Subject",
	      "NCBI_Description",
	      "NCBI_Bitscore",
	      "InterproScan result available",
	      "HMMER subject",
	      "HMMER bitscore",
	      "JACKHMMER subject",
	      "JACKHMMER bitscore"
    )."\n";

print FH $line;
$ctx_output->add($line);

foreach my $id (sort {$a <=> $b} (keys %mascot_ids))
{
    my $interpro_present = "NO";
    if ( @{$mascot_ids{$id}{interpro}} > 0 )
    {
	$interpro_present = sprintf("YES (%d result lines)", int(@{$mascot_ids{$id}{interpro}}));
	foreach (@{$mascot_ids{$id}{interpro}})
	{
	    $_ =~ s/^\S+_(\d+)(\.p\d+)/pitu_v1_$1$2/;
	    print INTERPRO $_;
	    $ctx_interpro->add($_);
	}
    }

    my @diffseq_fields = ();
    foreach my $diffseq_file (@deseq_results)
    {
	my ($experiment) = $diffseq_file =~ /([^_]+)\.mat$/;

	push(@diffseq_fields, map { $mascot_ids{$id}{diffseq}{$experiment}{$_} } (qw(log2foldchange significant padj)));
    }

    my $line = join("\t",
		  "pitu_v1_".$id,
		  $mascot_ids{$id}{quantification}{'./stringtie/3_8_dta.gtf_cov'},
		  $mascot_ids{$id}{quantification}{'./stringtie/3_8_dta.gtf_tpm'},
		  $mascot_ids{$id}{quantification}{'./stringtie/6_dta.gtf_cov'},
		  $mascot_ids{$id}{quantification}{'./stringtie/6_dta.gtf_tpm'},
		  $mascot_ids{$id}{quantification}{'./stringtie/7_dta.gtf_cov'},
		  $mascot_ids{$id}{quantification}{'./stringtie/7_dta.gtf_tpm'},
		  $mascot_ids{$id}{deseqnormalized}{venomgland},
		  $mascot_ids{$id}{deseqnormalized}{femalebody},
		  $mascot_ids{$id}{deseqnormalized}{malebody},
		  @diffseq_fields,
		  join(";", @{$mascot_ids{$id}{manual}}),
		  ( map { exists $mascot_ids{$id}{toxprot}{$_} ? $mascot_ids{$id}{toxprot}{$_} : "" } ("subject", "description", "bitscore")),
		  ( map { exists $mascot_ids{$id}{ncbi}{$_} ? $mascot_ids{$id}{ncbi}{$_} : "" } ("subject", "description", "bitscore")),
		  $interpro_present,
		  ((exists $mascot_ids{$id}{hmmer}{hmmerhit}) ? $mascot_ids{$id}{hmmer}{hmmerhit} : ""),
		  ((exists $mascot_ids{$id}{hmmer}{bitscore}) ? $mascot_ids{$id}{hmmer}{bitscore} : ""),
		  ((exists $mascot_ids{$id}{jackhmmer}{jackhmmerhit}) ? $mascot_ids{$id}{jackhmmer}{jackhmmerhit} : ""),
		  ((exists $mascot_ids{$id}{jackhmmer}{bitscore}) ? $mascot_ids{$id}{jackhmmer}{bitscore} : ""),
	)."\n";
    print FH $line;
    $ctx_output->add($line);
}
close(INTERPRO) || die "Unable to close interpro output table file '$interpro_out': $!\n";
close(FH) || die "Unable to close output table file '$table_out': $!\n";

$log->info(sprintf("Output file MD5: %s", $ctx_output->hexdigest));
$log->info(sprintf("Interpro output file MD5: %s", $ctx_interpro->hexdigest));
