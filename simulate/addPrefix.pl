#!/usr/bin/perl
use strict;
use warnings;

sub add(){
my ($file_1, $file_2, $k) = @_;
open IP,"< $file_1" or die $!;
open OP,"> $file_2" or die $!;
my $n = 0;
while(<IP>)
{
	
	chomp;
	if(/^\@FC:(\d+)/)	
	{
		my $str = "\@PRE_".$k."_";
		s/\@/$str/;
		if($n%4 != 0){ die "the quality line is recognized Qname!!!";}
	}
	print OP $_."\n";
	$n++;
}
close IP;
close OP;
}
&add($ARGV[0], $ARGV[1], $ARGV[2]);
