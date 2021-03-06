#!/usr/bin/env perl

use strict;
use warnings;

# input/output settings
my $compara = $ARGV[0] || 'Compara.homologies.34.tsv';
my $query = load_query($ARGV[1]);
my $species = load_species($ARGV[2]);
my $query_compara = $ARGV[3] || "YES";
my $prefix = $ARGV[4] || 'fungi_atg';

main();

# subroutines

sub main {
	if (uc($query_compara) eq 'YES') {
		query_protein_in_compara($compara, $query);
	}
	my $table = generate_ortholog_table($query, $species, $prefix);
	output_ortholog_table($table, $query, $prefix);
	return 1;
}

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
		$query{$w[0]}{$w[3]}{'gene'} = $w[1];
		$query{$w[0]}{$w[3]}{'protein'} = $w[2];
	}
	close IN;
	return \%query;
}

sub query_protein_in_compara {
	my ($compara, $query) = @_;
	system("mkdir -p $prefix");
	foreach my $id (sort by_string_number keys %{$query}) {
		unless (-e "$prefix/$id.txt") {
			foreach my $species (sort keys %{$query->{$id}}) {
				# caputre the lines matching query protein_id, probably we can use sql here later
				system("grep $query->{$id}->{$species}->{'protein'} $compara >> $prefix/$id.txt");
			}
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
		next if ($w[1] ne $protein);
		if (exists $ortholog{$w[7]}{'ortholog_gene_stable_id'}) {
			next if ($ortholog{$w[7]}{'ortholog_gene_stable_id'} =~ m/$w[5]/);
			$ortholog{$w[7]}{'ortholog_gene_stable_id'} .= ';' . $w[5];
		} else {
			$ortholog{$w[7]}{'ortholog_gene_stable_id'} = $w[5];
		}
		if (exists $ortholog{$w[7]}{'ortholog_protein_stable_id'}) {
			next if ($ortholog{$w[7]}{'ortholog_protein_stable_id'} =~ m/$w[6]/);
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

	my %table = ();
	foreach my $id (sort by_string_number keys %{$query}) {
		foreach my $species_name (sort keys %{$species}) {
			$table{$species_name}{$id}{'found'} = 'no';
			foreach my $query_species_name (sort keys %{$query->{$id}}) {
				my $protein = $query->{$id}->{$query_species_name}->{'protein'};
				my $ortholog = parse_homology_table("$prefix/$id.txt", $protein);
				if ($species_name eq $query_species_name) {
					next if ($table{$species_name}{$id}{'found'} eq 'yes');
					$table{$species_name}{$id}{'number'} = 1;
					$table{$species_name}{$id}{'ortholog_gene_id'} = $query->{$id}->{$query_species_name}->{'gene'};
					$table{$species_name}{$id}{'ortholog_protein_id'} = $query->{$id}->{$query_species_name}->{'protein'};
					$table{$species_name}{$id}{'found'} = 'yes';
				} elsif (exists $ortholog->{$species_name}) {
					if ($table{$species_name}{$id}{'found'} eq 'yes') {
						my ($ortholog_gene_id, $num1) = uniq_id($table{$species_name}{$id}{'ortholog_gene_id'}, $ortholog->{$species_name}->{'ortholog_gene_stable_id'});
						my ($ortholog_protein_id, $num2) = uniq_id($table{$species_name}{$id}{'ortholog_protein_id'}, $ortholog->{$species_name}->{'ortholog_protein_stable_id'});
						$table{$species_name}{$id}{'number'} = $num1;
						$table{$species_name}{$id}{'ortholog_gene_id'} = $ortholog_gene_id;
						$table{$species_name}{$id}{'ortholog_protein_id'} = $ortholog_protein_id;
					} else {
						$table{$species_name}{$id}{'number'} = $ortholog->{$species_name}->{'ortholog_number'};
						$table{$species_name}{$id}{'ortholog_gene_id'} = $ortholog->{$species_name}->{'ortholog_gene_stable_id'};
						$table{$species_name}{$id}{'ortholog_protein_id'} = $ortholog->{$species_name}->{'ortholog_protein_stable_id'};
						$table{$species_name}{$id}{'found'} = 'yes';
					}
				} else {
					next if ($table{$species_name}{$id}{'found'} eq 'yes');
					$table{$species_name}{$id}{'number'} = 0;
					$table{$species_name}{$id}{'ortholog_gene_id'} = '-';
					$table{$species_name}{$id}{'ortholog_protein_id'} = '-';
					$table{$species_name}{$id}{'found'} = 'no';
				}
			}
		}
	}
	return \%table;
}

sub output_ortholog_table {
	my ($table, $query, $prefix) = @_;

	open (NUMBER, ">>$prefix.number.tsv") or die "Cannot write $prefix.number.tsv : $!\n";
	open (TABLE, ">>$prefix.table.tsv") or die "Cannot write $prefix.table.tsv : $!\n";
	# print header
	foreach my $id (sort by_string_number keys %{$query}) {
		print NUMBER "\t", $id;
		print TABLE "\t", $id;
	}
	print NUMBER "\n";
	print TABLE "\n";
	# print body
	foreach my $species_name (sort keys %{$table}) {
		print NUMBER format_species_name($species_name);
		print TABLE format_species_name($species_name);
		foreach my $id (sort by_string_number keys %{$query}) {
			print NUMBER "\t", $table->{$species_name}->{$id}->{'number'};
			print TABLE "\t", $table->{$species_name}->{$id}->{'ortholog_gene_id'}, "/", $table->{$species_name}->{$id}->{'ortholog_protein_id'};
		}
		print NUMBER "\n";
		print TABLE "\n";
	}
	close TABLE;
	close NUMBER;
	return 1;
}

sub uniq_id {
	my ($str1, $str2) = @_;
	my @w1 = split /\;/, $str1;
	my @w2 = split /\;/, $str2;
	my %uniq = ();
	my $uniq = '';
	my $num = 0;
	foreach my $str (@w1) {
		$uniq{$str} = 1;
	}
	foreach my $str (@w2) {
		$uniq{$str} = 1;
	}
	foreach my $w (sort keys %uniq) {
		$uniq .= $w . ";";
		$num++;
	}
	$uniq =~ s/\;$//;
	return ($uniq, $num);
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
	print STDERR "$file exists, will remove in $time seconds\n";
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
	# specific case
	return "Dacryopinax sp. DJM-731 SS1" if $name eq 'dacryopinax_sp_djm_731_ss1';
	return "Nematocida sp. 1 ERTm2" if $name eq 'nematocida_sp_1_ertm2';
	return "Pyrenophora tritici-repentis" if $name eq 'pyrenophora_triticirepentis';
	return "Saccharomyces sp. 'boulardii'" if $name eq 'saccharomyces_sp_boulardii_';
	return "Saccharomycetaceae sp. 'Ashbya aceri'" if $name eq 'saccharomycetaceae_sp_ashbya_aceri_';
	# format species name
	my @w = split /\_/, $name;
	$w[0] =~ s/(\w+)/\u$1/;
	$name = $w[0] . ' ' . $w[1];
	return $name;
}
