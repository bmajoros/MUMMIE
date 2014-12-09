#!/usr/bin/perl
use strict;

my (@counts,$sum);
while(<STDIN>) {
  chomp;
  if(/get-extra/) {
    foreach my $pair (@counts) {
      my ($count,$type)=@$pair;
      my $percent=int(100*$count/$sum+5/9);
      print "$percent\%\t$type\n";
    }
    $sum=0;
    @counts=();
    undef @counts;
  }
  elsif(/(\d+)\s+(\S+)/) {
    my ($count,$type)=($1,$2);
    push @counts,[$count,$type];
    $sum+=$count;
  }
}

