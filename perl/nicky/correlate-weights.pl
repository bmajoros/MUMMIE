#!/usr/bin/perl
use strict;

my $file1="WT-weights.txt";
my $file2="D1-weights.txt";

my %weights1;
open(IN,$file1) || die;
while(<IN>) {
  if(/(\S+)\s(\S+)/) {
    $weights1{$1}=$2;
  }
}
close(IN);
open(IN,$file2) || die;
while(<IN>) {
  if(/(\S+)\s(\S+)/) {
    my $w1=$weights1{$1}+0;
    my $w2=$2;
    next unless $w1>0 && $w2>0;
    $w1=log($w1);
    $w2=log($w2);
    print "$w1 $w2\n";
  }
}
close(IN);

