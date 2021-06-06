#!/usr/bin/env perl

use warnings;
use strict;

# in order to search for modules in the directory where this script is located
use File::Basename;
use Cwd;
use lib dirname (&Cwd::abs_path(__FILE__));

# modules in this distribution
use GetOptions;

# names of command line options
my $HELP_OPTION = 'help';
my $FASTA_OPTION = 'fasta';
my $TAXONOMY_OPTION = 'taxonomy';
my $COLORS_OPTION = 'colors';
my $OUTPUT_OPTION = 'output';

# types for command line options; see 'Getopt::Long' Perl documentation for information on option types
my %OPTION_TYPES = ($HELP_OPTION => '',
		    $FASTA_OPTION => '=s',
		    $TAXONOMY_OPTION => '=s',
		    $COLORS_OPTION => '=s',
		    $OUTPUT_OPTION => '=s');

sub main 
{
    my %args = &GetOptions::get_options(\%OPTION_TYPES);

    if ($args{$HELP_OPTION} || (not ($args{$FASTA_OPTION} && $args{$TAXONOMY_OPTION} && $args{$COLORS_OPTION} &&
				     $args{$OUTPUT_OPTION})))
    {
	&help();
	exit;
    }

    open(TAXONOMY, $args{$TAXONOMY_OPTION}) || die "ERROR: could not open file '$args{$TAXONOMY_OPTION}' for read";
    open(COLORS, $args{$COLORS_OPTION}) || die "ERROR: could not open file '$args{$COLORS_OPTION}' for read";
    open(FASTA, $args{$FASTA_OPTION}) || die "ERROR: could not open file '$args{$FASTA_OPTION}' for read";
    open(OUTPUT, ">", $args{$OUTPUT_OPTION}) || die "ERROR: could not open file '$args{$OUTPUT_OPTION}' for write";

    my $taxonomy_header = <TAXONOMY>;

    chomp $taxonomy_header;
    
    my ($id_column, $name_column, $lineage_column) = (-1, -1, -1);

    my @taxonomy_header_columns = split("\t", $taxonomy_header);

    for (my $index = 0; $index < scalar(@taxonomy_header_columns); $index++)
    {
	if ($taxonomy_header_columns[$index] eq "EukProt_ID")
	{
	    $id_column = $index;
	}
	elsif ($taxonomy_header_columns[$index] eq "Name_to_Use")
	{
	    $name_column = $index;
	}
	elsif ($taxonomy_header_columns[$index] eq "Taxonomy_UniEuk")
	{
	    $lineage_column = $index;
	}
    }
    
    my %taxonomy_by_data_set = ();
    
    while (<TAXONOMY>)
    {
	chomp;

	my @taxonomy_columns = split("\t");
	
	$taxonomy_by_data_set{$taxonomy_columns[$id_column] . "_" . $taxonomy_columns[$name_column]} =
	    $taxonomy_columns[$lineage_column];
    }

    close TAXONOMY;

    my %colors_by_data_set = ();

    my $reading_colors = 0;
    
    while (<COLORS>)
    {
	if ($reading_colors)
	{
	    chomp;

	    my ($node, undef, $color, undef) = split(",");

	    foreach my $data_set (keys %taxonomy_by_data_set)
	    {
		if ($taxonomy_by_data_set{$data_set} =~ /$node/)
		{
		    $colors_by_data_set{$data_set} = join(" ", ($color, $node));
		}
	    }
	}
	elsif (/^DATA/)
	{
	    $reading_colors = 1;
	}
    }

    close COLORS;

    print OUTPUT "DATASET_COLORSTRIP\n";
    print OUTPUT "SEPARATOR SPACE\n";
    print OUTPUT "DATASET_LABEL Taxonomy\n";
    print OUTPUT "COLOR #ffffff\n";
    print OUTPUT "COLOR_BRANCHES 1\n";

    print OUTPUT "DATA\n";

    while (<FASTA>)
    {
	if (/^>/)
	{
	    chomp;

	    my ($record_name, undef) = split(/\s+/);

	    $record_name = substr($record_name, 1);

	    my ($data_set, undef) = split(/\|/, $record_name);

	    if (not defined $colors_by_data_set{$data_set})
	    {
		print $data_set . "\n";
	    }
	    
	    print OUTPUT join(" ", ($record_name, $colors_by_data_set{$data_set}, "mouseover")) . "\n";
	    
	}
    }
    
    close FASTA;
    close OUTPUT;
}

sub help
{
    my $HELP = <<HELP;
Syntax: $0 -$FASTA_OPTION <fasta> -$TAXONOMY_OPTION <txt> -$COLORS_OPTION <txt> -$OUTPUT_OPTION <txt>

Produce an annotation file for iTOL, of type "COLORSTRIP", which
colors sequences in a FASTA file by the taxonomic lineage of the
species in which they are found.

    -$HELP_OPTION : print this message
    -$FASTA_OPTION : FASTA file containing records to annotate
    -$TAXONOMY_OPTION : tab-delimited spreadsheet linking species to their taxonomic lineage (from EukProt)
    -$COLORS_OPTION : list of colors to be used (from EukProt)
    -$OUTPUT_OPTION : output file, to be used as input to iTOL

   
HELP

    print STDERR $HELP;
}

&main();
