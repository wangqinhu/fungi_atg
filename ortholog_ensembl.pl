#!/usr/bin/env perl

use strict;
use warnings;

# input/output settings
my $compara = $ARGV[0] || 'Compara.homologies.32.tsv';
my $query = load_query($ARGV[1]);
my $species = load_species($ARGV[2]);
my $prefix = $ARGV[3] || 'fungi_atg';

query_protein_in_compara($compara, $query);
generate_ortholog_table($query, $species, $prefix);

# subroutines
sub load_query {
	my $file = shift;
	my %query = ();
	open (IN, $file) or die "Cannot open $file : $!";
	while (<IN>) {
		chomp;
		next if /^\#/;
		next if /^\t/;
		next if /^\s*$/;
		my @w = split /\t/;
		$query{$w[0]}{'gene'} = $w[1];
		$query{$w[0]}{'protein'} = $w[2];
		$query{$w[0]}{'species'} = $w[3];
	}
	close IN;
	return \%query;
}

sub query_protein_in_compara {
	my ($compara, $query) = @_;
	system("mkdir -p $prefix");
	foreach my $id (sort by_string_number keys %{$query}) {
		unless (-e "$prefix/$id.txt") {
			system("grep $query->{$id}->{'gene'} $compara > $prefix/$id.txt");
		}
	}
	return 1;
}

sub load_species {
	my $file = shift;
	my %species = ();
	open (IN, $file) or die "Cannot open $file : $!\n";
	while (<IN>) {
		chomp;
		next if /^\s*$/;
		$species{$_} = 1;	
	}
	close IN;
	return \%species;
}

sub parse_homology_table {
	my ($file, $protein) = @_;
	my %ortholog = ();
	open (IN, $file) or die "Cannot open $file : $!\n";
	while (<IN>) {
		chomp;
		next if /paralog/;
		my @w = split /\t/;
		next if ($w[0] ne $protein);
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
	my ($query, $species, $prefix) = @_;

	if (-e "$prefix.table.tsv") {
		remove_file("$prefix.table.tsv", 5);
	}
	if (-e "$prefix.number.tsv") {
		remove_file("$prefix.number.tsv", 5);
	}

	foreach my $id (sort by_string_number keys %{$query}) {
		my $protein = $query->{$id}->{'protein'};
		my $ortholog = parse_homology_table("$prefix/$id.txt", $protein);
		open (TABLE, ">>$prefix.table.tsv") or die "Cannot write $prefix.table.tsv : $!\n";
		open (NUMBER, ">>$prefix.number.tsv") or die "Cannot write $prefix.number.tsv : $!\n";
		foreach my $species_name (sort keys %{$species}) {
			my $buffer_number = '';
			my $buffer_details = '';
			if (exists $ortholog->{$species_name}) {
				$buffer_number = $id . "\t" .
					format_species_name($species_name) . "\t" .
					$ortholog->{$species_name}->{'ortholog_number'};
				$buffer_details = 
					$ortholog->{$species_name}->{'homology_type'} . "\t" .
					$ortholog->{$species_name}->{'ortholog_gene_stable_id'} . "\t" .	
					$ortholog->{$species_name}->{'ortholog_protein_stable_id'};
			} else {
				if ($species_name eq $query->{$id}->{'species'}) {
					$buffer_number = $id . "\t" . format_species_name($species_name) . "\t" . "1";
					$buffer_details = "query" . "\t" . $query->{$id}->{'gene'} . "\t" . $query->{$id}->{'protein'};
				} else {
					$buffer_number = $id . "\t" . format_species_name($species_name) . "\t" . "0";
					$buffer_details = "na" . "\t-\t-";
				}
			}
			print NUMBER $buffer_number, "\n";
			print TABLE $buffer_number, "\t", $buffer_details, "\n";
		}
		close TABLE;
		close NUMBER;
	}
	return 1;
}

sub by_string_number {
	$a =~ /(\d+)/;
	my $numa = $1;
	$b =~ /(\d+)/;
	my $numb = $1;
	return $numa <=> $numb;
}

sub remove_file {
	my ($file, $time) = @_;
	$time = 5 if (!defined $time);
	print STDERR "$file exists, will remove in $time seconds ... \n";
	for my $i (1..$time) {
		print STDERR ".";
		sleep(1);
	}
	system("rm $file");
	print STDERR "\nRemoved.\n";
	return 1;
}

sub format_species_name {
	my $name =  shift;
	my @w = split /\_/, $name;
	$w[0] =~ s/(\w+)/\u$1/;
	$name = $w[0] . ' ' . $w[1];
	return $name;
}
