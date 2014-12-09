#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;

my $SCUMM="/home/ohler/bmajoros/SCUMM";
my $BASE="$SCUMM/test/rand";

my @sums;

my @files=`ls $BASE/*-trained.hmm`;
foreach my $file (@files) {
  chomp $file;
  my $acc=`$SCUMM/test.pl $file test`;
  $acc=~/(.*)\%/ || die;
  $acc=$1;
  $file=~/rand(\d+)-trained/ || die $file;
  my $id=$1;
  open(IN,"$BASE/rand$id.stdout") || die;
  my ($initialLL,$finalLL);
  while(<IN>) {
    if(/ITERATION #(\d+)\s+LL=(-\d+)/) {
      my ($iter,$LL)=($1,$2);
      if($iter==0) {$initialLL=$LL}
      elsif($iter==99) {$finalLL=$LL}
      $sums[$iter]+=$LL;
    }
  }
  close(IN);
  print "$id\t$initialLL\t$finalLL\t$acc\n";
}

open(OUT,">convergence.txt") || die;
for(my $i=0 ; $i<100 ; ++$i) {
  my $ave=$sums[$i]/100;
  my $iter=$i+1;
  print OUT "$iter $ave\n";
}
close(OUT);



