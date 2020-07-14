#!/usr/bin/perl
use strict;
use warnings;

my $output = $ARGV[-1];
my $inputnum = scalar(@ARGV) -1;
print $inputnum."\n";
my $str = "cat ";
for my $i (0..($inputnum -1))
{
    $str .= " ".$ARGV[$i];
}
$str .= "> $output";
#print $str;
system  $str;
