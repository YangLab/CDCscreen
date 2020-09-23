#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

sub usage {
	print <<"END_USAGE";
Usage: perl $0
Needed:
	--input input file,sep with ","
	--output output file
	--suffix example:txt or vcf or ...
Optional:
	--col default ALL!!! Start col=0,sep with "/," exp 1/2/3,3/4,6/7,1
	--head if file has header exp 0,1,0,1
	--fill default 0
END_USAGE
	exit;
}

my ($input,$output,$suffix,$head,$fill,$col);
GetOptions (
	'input=s'=>\$input,
	'output=s'=>\$output,
	'suffix=s'=>\$suffix,
	'col=s'=>\$col,
	'fill=s'=>\$fill,
	'head=s'=>\$head,
) or usage();
usage() if (!$input or !$output or !$suffix);

my @file=split(/,/,$input);
if (defined $col){
	my @col=split(/,/,$col);
	die "sample counts!= col counts\n" if ($#col!=$#file);
}
if (defined $head){
	my @headcount=split(/,/,$head);
	die "sample counts!= head counts\n "if ($#headcount!=$#file);
}
$fill=0 if (!$fill);
my %hash;
my @sample;
for my $num(0..$#file){
	open my $file_in,"$file[$num]";
	$file[$num]=~s/\.$suffix$// if ($suffix);
	if (defined $head){
		my @headcount=split(/,/,$head);
		<$file_in> if ($headcount[$num]==1);
	}
	while (<$file_in>){
		chomp;
		my @F=split;
		my $value;
		if (defined $col){
			my @col=split(/,/,$col);
			my @col_sample=split(/\//,$col[$num]);
			if ($#col_sample==0){
				$value=$F[$col_sample[0]];
			}
			elsif($#col_sample>=0){
				$value=$F[$col_sample[0]];
				$value=join "\t",$value,$F[$col_sample[$_]] for (1..$#col_sample);
			}
		}
		else{
			$value=join "\t",@F[1..$#F];
		}
		$hash{$F[0]}{$file[$num]}=$value;
	}
	close $file_in;
	push @sample,$file[$num];
}
open my $output_file,">$output";
print $output_file "Event";
print $output_file "\t$_" for @sample;
print $output_file "\n";
for my $keys (keys %hash){
	print $output_file "$keys";
	for my $sample (@sample){
		print $output_file "\t$hash{$keys}{$sample}" if (exists $hash{$keys}{$sample});
		print $output_file "\t$fill" if (!exists $hash{$keys}{$sample});
	}
	print $output_file "\n";
}
close $output_file;
