#!/usr/bin/perl
use strict;
use warnings;

if(scalar(@ARGV)%2 ne 0) {die "error num of arguments";}
my $num = scalar(@ARGV/2);
my @vcf = ();
my @proportion = ();
for my $n(0..($num-1)){$vcf[$n] = $ARGV[$n]; $proportion[$n] = $ARGV[$n+$num];}
#print "@vcf"."\n";
#print "@proportion"."\n";


sub sumOfBase(){
	my ($tabfile, $pro) = @_;
	open IP,"< $tabfile" or die $!;
	my $length = 0;
	while(<IP>)
	{
		chomp;
		next if(/^\#/);
		my ($start , $end, $copy_num) = (split/\t/,$_)[1,2,3];
		$length += ($copy_num - 1) * ($end - $start);
	}
	$length = $length * $pro;
	close IP;
	return $length;
}
my $totalLength = 0;
for my $n(0..($num-1)){ $totalLength += &sumOfBase($vcf[$n], $proportion[$n]);}
open OP,"> ploidy.txt" or die $!;
#whole length : 3137454505
#length chr1 - chr22 : 2881033286 
my $total = 2881033286;
print OP "CNA_length\thaploid_length\tploidy\n";
print OP $totalLength."\t".$total."\t".($totalLength/$total + 2)."\n";
close OP;
