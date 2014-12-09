#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;

my %names;
open(IN,"iulian-names.txt") || die;
while(<IN>) {
  chomp;
  if(/(\S+)\s+(\S+)/){$names{$1}=$2}
}
close(IN);

my @clusters;
my @files=`ls iulian/*.txt`;
foreach my $file (@files) {
  chomp $file;
  my $cluster=[];
  open(IN,$file) || die;
  while(<IN>) {
    if(/(\d+)/) {push @$cluster,$names{$1}}
  }
  close(IN);
  push @clusters,$cluster;
}

@clusters=sort {0+@$b <=>  0+@$a} @clusters;
my $n=@clusters;
for(my $i=0 ; $i<$n; ++$i) {
  my $cluster=$clusters[$i];
  print "cluster $i: ";
  foreach my $elem (@$cluster) {print "$elem\t"}
  print "\n";
}
