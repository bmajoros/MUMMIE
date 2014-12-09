#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;

my %hash;
open(IN,"UTRs.txt") || die;
while(<IN>) {
  chomp;
  my @fields=split/\s+/,$_;
  my ($chr,$begin,$end,$gene)=@fields;
  push @{$hash{$chr}},[$begin,$end];
}
close(IN);

my @keys=keys %hash;
foreach my $key (@keys) {
  my $array=$hash{$key};
  my @sorted=sort {$a->[0] <=> $b->[0]} @$array;
  my $n=@sorted;
  for(my $i=0 ; $i<$n-1 ; ++$i) {
    my $this=$sorted[$i];
    my $next=$sorted[$i+1];
    my ($thisBegin,$thisEnd)=@$this;
    my ($nextBegin,$nextEnd)=@$next;
    if($thisBegin<$nextEnd && $nextBegin<$thisEnd) {
      die "OVERLAP\n$thisBegin $thisEnd $nextBegin $nextEnd\n"}
  }
}


