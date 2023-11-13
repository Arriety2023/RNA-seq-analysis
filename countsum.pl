#!/usr/bin/perl
use strict;use warnings;
my ($num,$i,$n);
$num=0;
$i=0;
$n=0;
while (<>) {
$i++;
$n++;
if ($i eq 2) {
	chomp($_);
	$num=$num+length($_);
	$i=-2;	
}else {
next;  
}
}
$n=$n/4;
print "$n\t$num";
