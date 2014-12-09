#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;

my @A;
open(IN,"tmp.99") || die;
while(<IN>) {
  chomp;
  if(/(\S+)\s+(\S+)/){push(@A,[$1,$2])}
}
close(IN);
die unless @A==4;

open(IN,"tmp.100") || die;
while(<IN>) {
  chomp;
  if(/(\S.*\S)/) {
    my @fields=split/\s+/,$_;
    next unless @fields==4;
    my ($a,$b);
    for(my $i=0 ; $i<4 ; ++$i) {
      $a+=$fields[$i]*$A[$i]->[0];
      $b+=$fields[$i]*$A[$i]->[1];
    }
    print "$a\t$b\n";
  }
}
close(IN);

