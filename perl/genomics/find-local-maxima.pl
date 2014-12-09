#!/usr/bin/perl
use strict;

my $usage="$0 <file>";
die "$usage\n" unless @ARGV==1;
my ($filename)=@ARGV;

my ($prevY,$prevSlope);
open(IN,$filename) || die "can't open $filename\n";
while(<IN>)
  {
    if(/(\S+)\s+(\S+)/)
      {
	my ($x,$y)=($1,$2);
	if(!defined($prevY)) {$prevY=$y}
	my $slope=$y-$prevY;
	if($prevSlope>0 && $slope<0)
	  {
	    print "$x $y\n";
	  }
	$prevSlope=$slope;
	$prevY=$y;
      }
  }
close(IN);




