#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <scores.txt> <cutoff>\n" unless @ARGV==2;
my ($infile,$cutoff)=@ARGV;

open(IN,$infile) || die "can't open $infile";
my ($above,$below);
while(<IN>) {
  chomp;
  my $x=$_+0;
  if($x<$cutoff) {++$below}
  else {++$above}
}
close(IN);
my $FNrate=$below/($below+$above);
print "$FNrate\n";

