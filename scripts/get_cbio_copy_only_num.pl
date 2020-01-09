#! usr/bin/env perl
use strict;
use warnings FATAL => 'all';

my $in=shift;
my $in2=shift;
my %hash;
my %cbio;

open IN, $in;
while(<IN>){
	chomp;
	my @a=split(/\t/,$_);
	$hash{$a[0]}=$a[1];
	}
open IN2, $in2;
while(<IN2>){
	chomp;
	my @b=split(/\t/,$_);
	my $num=$b[3];
	if (exists $hash{$b[-1]}){
		my $as="$b[-1]\t$hash{$b[-1]}";
#			print "$b[-1]\t$hash{$b[-1]}\t$num\n";
			$cbio{$as}=$num;
			}
	}
foreach my $key (sort keys %cbio){
	print "$key\t$cbio{$key}\n";
		}
