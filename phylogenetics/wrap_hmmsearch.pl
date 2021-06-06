#!/usr/bin/env perl

use warnings;
use strict;

# core Perl modules
use IPC::Run;

# in order to search for modules in the directory where this script is located
use File::Basename;
use Cwd;
use lib dirname (&Cwd::abs_path(__FILE__));

# modules in this distribution
use GetOptions;

my $HELP_OPTION = 'help';
my $TARGET_OPTION = 'target';
my $OUTPUT_OPTION = 'output';
my $HMM_OPTION = 'hmm';
my $EVALUE_OPTION = 'evalue';
my $HITS_OPTION = 'hits';

my %OPTION_TYPES = ($HELP_OPTION => '',
		    $TARGET_OPTION => '=s',
		    $OUTPUT_OPTION => '=s',
		    $HMM_OPTION => '=s',
		    $EVALUE_OPTION => '=f',
		    $HITS_OPTION => '=i');

my $DIRECTORY_DELIMITER = "/";
my $TEMPORARY_HMMER_OUTPUT = "/tmp/hmmer.out";
my $DEFAULT_HITS = 1;

sub main 
{
    my %args = &GetOptions::get_options(\%OPTION_TYPES);

    if ($args{$HELP_OPTION} || (not ($args{$HMM_OPTION} && $args{$OUTPUT_OPTION} && $args{$TARGET_OPTION})))
    {
	&help();
	exit;
    }

    my $hits_to_retrieve = $args{$HITS_OPTION} || $DEFAULT_HITS;
    
    my %hmms_to_search = ();

    if (-d $args{$HMM_OPTION})
    {
	if (not -d $args{$OUTPUT_OPTION})
	{
	    die "if a directory is specified in '-$HMM_OPTION', then '-$OUTPUT_OPTION' must also be a directory, and " .
		"'$args{$OUTPUT_OPTION}' is not a directory";
	}

	opendir(HMM_DIR, $args{$HMM_OPTION}) || die "could not open directory '$args{$HMM_OPTION}' for read";

	foreach my $file (readdir(HMM_DIR))
	{
	    if ($file =~ /(.*)\.hmm$/)
	    {
		my $hmm_name = $1;

		my $hmm_file_path = $args{$HMM_OPTION} . $DIRECTORY_DELIMITER . $file;

		my $output_file_path = $args{$OUTPUT_OPTION} . $DIRECTORY_DELIMITER . $hmm_name . ".fasta";

		open($hmms_to_search{$hmm_file_path}, ">", $output_file_path) || die "could not open " .
		    "'$output_file_path' for write";
	    }
	}

	closedir HMM_DIR;
    }
    else
    {
	open($hmms_to_search{$args{$HMM_OPTION}}, ">", $args{$OUTPUT_OPTION}) || die "could not open " .
	    "'$args{$OUTPUT_OPTION}' for write";
    }
		
    my %fasta_files = ();

    if (-d $args{$TARGET_OPTION})
    {
	opendir(DIR, $args{$TARGET_OPTION}) || die "could not open directory '$args{$TARGET_OPTION}' for read";
	
	foreach my $file (readdir(DIR))
	{
	    if ($file =~ /(.*?)\.fasta$/)
	    {
		$fasta_files{$1} = $args{$TARGET_OPTION} . $DIRECTORY_DELIMITER . $file;
	    }
	}

	closedir DIR;
    }
    else
    {
	if ($args{$TARGET_OPTION} =~ /(.*?)\.fasta$/)
	{
	    die "ERROR: could not parse FASTA file name '$args{$TARGET_OPTION}'";
	}

	$fasta_files{$1} = $args{$TARGET_OPTION};
    }

    foreach my $species (sort keys %fasta_files)
    {
	print STDOUT "--- searching $species ---\n";

	# get the sequences of all of the target proteins
	my %target_sequences = ();
	
	open(INPUT, $fasta_files{$species}) || die "could not open '$fasta_files{$species}' for read";
	
	my $current_id = "";

	while (<INPUT>)
	{
	    chomp;
	    
	    if (/^>(.*)/)
	    {
		$current_id = join(" ", (split(/\s+/, $1)));
	    }
	    else
	    {
		$target_sequences{$current_id} .= $_;
	    }
	}

	close INPUT || die "could not close '$fasta_files{$species}' after read";

	foreach my $hmm_file (keys %hmms_to_search)
	{
	    my %gene_ids_already_printed = ();
	    
	    my ($ids, $e_values, $starts, $ends, $descriptions) =
		&run_and_parse_hmmer($hmm_file, $fasta_files{$species}, $args{$EVALUE_OPTION}, $hits_to_retrieve);
	    
	    if (defined $ids->[0])
	    {
		my %domain_counter = ();
		
		for (my $index = 0; $index < scalar(@$ids); $index++)
		{
		    my $hit_id = join(" ", ($ids->[$index], $descriptions->[$index]));
		    
		    if ($descriptions->[$index] eq "-")
		    {
			$hit_id = $ids->[$index];
		    }
		    
		    my $target_sequence = $target_sequences{$hit_id};
		    
		    if (not defined $target_sequence)
		    {
			die "no sequence for '$hit_id' in target database '$fasta_files{$species}'";
		    }
		    
		    my $gene_id = $ids->[$index];

		    $gene_id =~ s/[:;]/_/g;
		
		    if (not defined $gene_ids_already_printed{$gene_id})
		    {
			$gene_ids_already_printed{$gene_id} = 1;
		    }
		    else
		    {
			next;
		    }
		    
		    if (scalar(keys %fasta_files) > 1)
		    {
#			print { $hmms_to_search{$hmm_file} } ">" . $species . "|" . $gene_id . "_" . 
#			    ($index + 1) . "_" . $e_values->[$index] . "\n";

			print { $hmms_to_search{$hmm_file} } ">" . $species . "|" . $gene_id . "\n";
			
			print { $hmms_to_search{$hmm_file} } $target_sequence . "\n";
		    }
		    else
		    {
			if (not defined $domain_counter{$gene_id})
			{
			    $domain_counter{$gene_id} = 0;
			}
			
			$domain_counter{$gene_id}++;
			
			print { $hmms_to_search{$hmm_file} } ">" . $gene_id . " " . $domain_counter{$gene_id} . " " .
			    $e_values->[$index] . "\n";
			
			print { $hmms_to_search{$hmm_file} } $target_sequence . "\n";
		    }
		}
	    }

	    undef $ids; undef $e_values; undef $starts; undef $ends;
	    undef %gene_ids_already_printed;
	}
	
	undef %target_sequences;
    }

    foreach my $hmm_file_path (keys %hmms_to_search)
    {
	close($hmms_to_search{$hmm_file_path}) || die "could not close output file for '$hmm_file_path' after write";
    }
}

sub run_and_parse_hmmer
{
    my ($hmm_file, $database, $e_value, $hits) = @_;

    if (not -e $hmm_file)
    {
	die "HMM '$hmm_file' does not exist";
    }

    if (not -e $database)
    {
	die "database '$database' does not exist";
    }

    my $output_format = "--tblout";

    my $e_value_threshold = "--cut_ga";
    
    if ($e_value)
    {
	$e_value_threshold = "-E " . $e_value;
    }

    my @hmmsearch_command = ("hmmsearch",
			     $e_value_threshold,
			     $output_format,
			     $TEMPORARY_HMMER_OUTPUT,
			     $hmm_file,
			     $database);

    my $hmmsearch_command = join(" ", @hmmsearch_command); 
    
    print STDOUT "$hmmsearch_command\n";
    
    &IPC::Run::run(\@hmmsearch_command, ">", "/dev/null") || die "ERROR: could not run '$hmmsearch_command': $?";

    open(INPUT, $TEMPORARY_HMMER_OUTPUT) || die "could not open '$TEMPORARY_HMMER_OUTPUT' for read";
    
    my @ids = ();
    my @e_values = ();
    my @starts = ();
    my @ends = ();
    my @descriptions = ();

    while (<INPUT>)
    {
	if (/^#/)
	{
	    next;
	}
	else
	{
	    chomp;

	    my @result_columns = split(/\s+/);
	    
	    push(@ids, $result_columns[0]);
	    push(@e_values, $result_columns[4]);
	    
	    my @output_descriptions = ();
	    
	    for (my $index = 18; $index < scalar(@result_columns); $index++)
	    {
		push(@output_descriptions, $result_columns[$index]);
	    }
	    
	    push(@descriptions, join(" ", @output_descriptions));
	    
	    undef @result_columns;
	    undef @output_descriptions;
	}
	
	if ($hits && scalar(@ids) == $hits)
	{
	    last;
	}
    }
    
    close INPUT;
	    
    return \@ids, \@e_values, \@starts, \@ends, \@descriptions;
}

sub help
{
    my $HELP = <<HELP;
Syntax: $0 -$HMM_OPTION <hmm or dir> -$OUTPUT_OPTION <FASTA or dir> -$TARGET_OPTION <FASTA or dir>
        [-$HITS_OPTION <int>] [-$EVALUE_OPTION <float>]
	
Search the input HMM (or directory of files ending in .hmm) against a
FASTA file (or a directory of files ending in .fasta) and retrieve the
sequence of each of the domain hits.

Reports i-Evalue and uses envelope coordinates of domain hits.

    -$HELP_OPTION : print this message
    -$HMM_OPTION : file or directory of HMM(s) to search with
    -$OUTPUT_OPTION : output FASTA file (or directory, if a directory is specified in '-$HMM_OPTION')
    -$TARGET_OPTION : target protein FASTA file or directory of FASTA files
    -$HITS_OPTION : return this many hits (starting with the best hit; default $DEFAULT_HITS)
    -$EVALUE_OPTION : instead of using Pfam gathering threshold, use this E value
    
HELP

    print STDERR $HELP;
}

&main();
