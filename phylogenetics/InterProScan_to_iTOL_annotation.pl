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
my $INTERPRO_OPTION = 'interpro';
my $OUTPUT_OPTION = 'output';

# types for command line options; see 'Getopt::Long' Perl documentation for information on option types
my %OPTION_TYPES = ($HELP_OPTION => '',
		    $INTERPRO_OPTION => '=s',
		    $OUTPUT_OPTION => '=s');

my @SHAPES = qw(RE HH HV EL DI TR TL PL PR PU PD OC GP);
# color palette from ColorBrewer <https://colorbrewer2.org>
my @COLORS = qw(a6cee3 1f78b4 b2df8a 33a02c fb9a99 e31a1c fdbf6f ff7f00 cab2d6 6a3d9a ffff99 b15928);

sub main 
{
    my %args = &GetOptions::get_options(\%OPTION_TYPES);

    if ($args{$HELP_OPTION} || (not ($args{$INTERPRO_OPTION} && $args{$OUTPUT_OPTION})))
    {
	&help();
	exit;
    }

    open(INTERPRO, $args{$INTERPRO_OPTION}) || die "ERROR: could not open file '$args{$INTERPRO_OPTION}' for read";
    open(OUTPUT, ">", $args{$OUTPUT_OPTION}) || die "ERROR: could not open file '$args{$OUTPUT_OPTION}' for write";

    my %domain_shapes = ();
    my %domain_colors = ();

    my $shape_index = 0;
    my $color_index = 0;
    
    my $legend_shapes = "LEGEND_SHAPES";
    my $legend_colors = "LEGEND_COLORS";
    my $legend_labels = "LEGEND_LABELS";

    my %domains_by_node = ();
    
    while (<INTERPRO>)
    {
	my ($node, undef, $length, undef, undef, $domain, $start, $end, undef) = split("\t");

	if (not defined $domain_shapes{$domain})
	{
	    $domain_shapes{$domain} = $SHAPES[$shape_index];
	    $domain_colors{$domain} = "#" . $COLORS[$color_index];

	    $legend_shapes .= "\t" . $SHAPES[$shape_index];
	    $legend_colors .= "\t" . "#" . $COLORS[$color_index];
	    $legend_labels .= "\t" . $domain;
	    
	    $shape_index++;

	    if ($shape_index == scalar(@SHAPES))
	    {
		$shape_index = 0;
	    }
	    
	    $color_index++;
	    
	    if ($color_index == scalar(@COLORS))
	    {
		$color_index = 0;
	    }
	}

	if ($node =~ /^(.*?)[:,]/)
	{
	    $node = $1;
	}

	if (not defined $domains_by_node{$node})
	{
	    $domains_by_node{$node} = [$length];
	}

	push(@{$domains_by_node{$node}}, join("|", ($domain_shapes{$domain}, $start, $end, $domain_colors{$domain}, $domain)));
    }

    print OUTPUT "DATASET_DOMAINS\n";
    print OUTPUT "SEPARATOR TAB\n";
    print OUTPUT "DATASET_LABEL\tDomains\n";
    print OUTPUT "COLOR\t#ffffff\n";

    print OUTPUT "LEGEND_TITLE\tDomains\n";
    print OUTPUT $legend_shapes . "\n";
    print OUTPUT $legend_colors . "\n";
    print OUTPUT $legend_labels . "\n";

    print OUTPUT "MARGIN\t25\n";
    print OUTPUT "SHOW_DOMAIN_LABELS\t0\n";
    
    print OUTPUT "DATA\n";
    
    foreach my $node (sort keys %domains_by_node)
    {
	print OUTPUT join("\t", ($node, @{$domains_by_node{$node}})) . "\n";
    }
    
    close INTERPRO;
    close OUTPUT;
}

sub help
{
    my $HELP = <<HELP;
Syntax: $0 -$INTERPRO_OPTION <tsv> -$OUTPUT_OPTION <txt>

Convert the output of InterProScan into annotations for use with iTOL.

    -$HELP_OPTION : print this message
    -$INTERPRO_OPTION : tab-delimited output of InterProScan (file ending in .tsv)
    -$OUTPUT_OPTION : output file, to be used as input to iTOL

   
HELP

    print STDERR $HELP;
}

&main();
