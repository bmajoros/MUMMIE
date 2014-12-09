#!/usr/bin/perl

my $offset=0;
if(@ARGV>0) {$offset=0+shift @ARGV}

my $x=0;
while(<STDIN>)
  {
    if(/\d/)
      {
	my $xc=$x+$offset;
	print "$xc $_";
	++$x;
      }
  }

