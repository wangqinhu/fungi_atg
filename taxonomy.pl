#!/usr/bin/perl

use strict;
use warnings;
use Bio::DB::Taxonomy;
use Bio::Tree::Tree;
use Bio::TreeIO;

# fungi
my @tax = load_name($ARGV[0]);

# taxonomy tree
tax2tree(@tax);
taxonomy(@tax);

# based on bioperl script bp_taxonomy2tree.pl 
sub tax2tree {
	my @species = @_;

	my $nodesfile = "data/taxonomy/nodes.dmp";
	my $namesfile = "data/taxonomy/names.dmp";

	my $db = Bio::DB::Taxonomy->new(-source => 'flatfile',
		-nodesfile => $nodesfile,
		-namesfile => $namesfile);

	# the full lineages of the species are merged into a single tree
	my $tree = undef;
	for my $name (@species) {
		my $ncbi_id = $db->get_taxonid($name);
		if ($ncbi_id) {
			my $node = $db->get_taxon(-taxonid => $ncbi_id);
			if ($tree) {
				$tree->merge_lineage($node);
			} else {
				$tree = new Bio::Tree::Tree(-node => $node);
			}
		} else {
			warn "no NCBI Taxonomy node for species ",$name,"\n";
		}
	}

	# simple paths are contracted by removing degree one nodes
	$tree->contract_linear_paths;

	# convert tree ids to their names for nice output with TreeIO
	foreach my $node ($tree->get_nodes) {
		$node->id($node->node_name);
	}

	# the tree is output in Newick format
	my $output = new Bio::TreeIO(-format => 'newick', -file => ">$ARGV[0].tree");
	$output->write_tree($tree);
}

# extract taxonomy information
sub taxonomy {
	my @taxnames = @_;
	my %lineages = ();

	my $nodesfile = "data/taxonomy/nodes.dmp";
	my $namesfile = "data/taxonomy/names.dmp";

	my $db = Bio::DB::Taxonomy->new(-source => 'flatfile',
		-nodesfile => $nodesfile,
		-namesfile => $namesfile);

	# based on:
	# https://github.com/gawbul/bioinformatics-scripts/blob/master/tax_identifier.pl
	foreach my $taxname (@taxnames) {
		my $taxon = $db->get_taxon(-name => $taxname);
		if (defined($taxon)) {
			my $tree = Bio::Tree::Tree->new(-node => $taxon);
			my @taxa = $tree->get_nodes;
			my $nodes = "";
			for (my $i = 1; $i < @taxa - 1; $i++) {
				my $t = $taxa[$i];
				my $nname = $db->get_taxon(-taxonid => $t->id());
				$nodes .= $nname->scientific_name() . "|";
			}
			$lineages{$taxname} = $nodes . $taxon->scientific_name();
		}
	}

	open (OUT, ">$ARGV[0].taxonomy.txt") or die "Cannot open file: $!\n";
	foreach my $taxname (sort @taxnames) {
		print OUT $taxname, "\t", $lineages{$taxname}, "\n";
	}
	close OUT;
}

sub load_name {
	my $file = shift;
	my @names = ();
	open (IN, $file) or die "Cannot open $file: \n";
	while (<IN>) {
		chomp;
		next if /^\#/;
		next if /^\s*$/;
		push @names, $_;
	}
	return @names;
}
