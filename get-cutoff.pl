#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <scores.txt> <FP-rate>\n" unless @ARGV==2;
my ($infile,$fp)=@ARGV;

open(IN,$infile) || die "can't open $infile";
my @array;
while(<IN>) {
  chomp;
  push @array,$_;
}
close(IN);
@array=sort {$a <=> $b} @array;
my $n=@array;
my $fpN=$fp*$n;
my $val=$array[$n-$fpN];
print "$val\n";

