#!/usr/bin/perl -w

die "please input two file " if @ARGV !=2;

open A,"$ARGV[0]" or die "$ARGV[0] can't be opened";
open B,"$ARGV[1]" or die "$ARGV[1] can't be opened";

my %INFO;

<A>;
while (<A>) {
	chomp;
	my @array=split "\t",$_;
	my $ID=$array[0];
	splice(\@array,0,1);
	my $p=join "\t",@array;
	$INFO{$ID}=$p;
}
while (<B>) {
	chomp;
	my @array=split "\t",$_;
	if (not exists $INFO{$array[0]} ) {$INFO{$array[0]}="No annotation"};
	print "$_\t$INFO{$array[0]}\n";
}
