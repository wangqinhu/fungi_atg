#!/usr/bin/env perl

use strict;
use warnings;

my $gene_id = $ARGV[0];
my $species = load_species($ARGV[1]);
my $ortholog = parse_homology_table($ARGV[2]);
generate_ortholog_table($gene_id, $species, $ortholog, $ARGV[3]);

# subroutines
sub load_species {
	my $file = shift;
	my %species = ();
	open (IN, $file) or die "Cannot open $file: $!\n";
	while (<IN>) {
		chomp;
		next if /^\s*$/;
		$species{$_} = 1;	
	}
	close IN;
	return \%species;
}

sub parse_homology_table {
	my $file = shift;
	my %ortholog = ();
	open (IN, $file) or die "Cannot open $file: $!\n";
	while (<IN>) {
		chomp;
		next if /paralog/;
		my @w = split /\t/;
		next if ($w[0] ne $gene_id);
		$ortholog{$w[7]}{'protein_stable_id'} = $w[1];
		$ortholog{$w[7]}{'homology_type'} = $w[3];
		if (exists $ortholog{$w[7]}{'ortholog_gene_stable_id'}) {
			$ortholog{$w[7]}{'ortholog_gene_stable_id'} .= ';' . $w[5];
		} else {
			$ortholog{$w[7]}{'ortholog_gene_stable_id'} = $w[5];
		}
		if (exists $ortholog{$w[7]}{'ortholog_protein_stable_id'}) {
			$ortholog{$w[7]}{'ortholog_protein_stable_id'} .= ';' . $w[6];
		} else {
			$ortholog{$w[7]}{'ortholog_protein_stable_id'} = $w[6];
		}
		$ortholog{$w[7]}{'ortholog_number'}++;
	}
	close IN;
	return \%ortholog;
}

sub generate_ortholog_table {
	my ($gene_id, $species, $ortholog, $file) = @_;
	open (OUT, ">$file") or die "Cannot write $file: $!\n";
	foreach my $species_name (sort keys %{$species}) {
		if (exists $ortholog->{$species_name}->{'homology_type'}) {
			my $buffer = $gene_id . "\t" .
			$ortholog->{$species_name}->{'protein_stable_id'} . "\t" .
			format_species_name($species_name) . "\t" .
			$ortholog->{$species_name}->{'ortholog_number'} . "\t" .
		   	$ortholog->{$species_name}->{'homology_type'} . "\t" .
			$ortholog->{$species_name}->{'ortholog_gene_stable_id'} . "\t" .	
			$ortholog->{$species_name}->{'ortholog_protein_stable_id'} . "\n";
			print OUT $buffer;
		} else {
			my $buffer = $gene_id . "\t-\t" . format_species_name($species_name) . "\t0\t-\t-\t-\n";
			print OUT $buffer;
		}
	}
	close OUT;
}

sub format_species_name {
	my $name =  shift;
	my @w = split /\_/, $name;
	$w[0] =~ s/(\w+)/\u$1/;
	$name = $w[0] . ' ' . $w[1];
	return $name;
}
