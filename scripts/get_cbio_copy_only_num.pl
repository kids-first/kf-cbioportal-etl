#! /usr/bin/env perl
# Script to output a headerless TSV with hugo gene symbol <tab> entrez ID <tab> copy number for downstrema processing
# Positional args, file 1 is a Hugo <tab> entrez ID TSV, file 2 is output from cnv_1_genome2gene.py
use strict;
use warnings FATAL => 'all';
use Try::Tiny;

my $in=shift;
my $in2=shift;
my %hash;
my %cbio;

# First file is hugo tsv, with gene symbols and entrez IDs
open IN, $in;
while(<IN>){
	chomp;
	my @a=split(/\t/,$_);
	if (@a == 2){
		$hash{$a[0]}=$a[1];
	}
	# Some genes have no Entrez ID, set to 0
	else{
		$hash{$a[0]} = 0;
	}
}
# Second file is output of bedtools intersect of ControlFreeC p value file filtered by cnv_1_genome2gene.py and a bed file with coords and gene names
open IN2, $in2;
while(<IN2>){
	try{
		chomp;
		my @b=split(/\t/,$_);
		# CN value
		my $num=$b[3];
		# ensure gene name is in hugo file
		if (exists $hash{$b[-1]}){
			# hashes by gene name <tab> entrez id 
			my $as="$b[-1]\t$hash{$b[-1]}";
				$cbio{$as}=$num;
				}
		}
		catch{
			warn "error thrown: $_ for file $in2\n";
			die;
		}
	}
foreach my $key (sort keys %cbio){
	print "$key\t$cbio{$key}\n";
		}
