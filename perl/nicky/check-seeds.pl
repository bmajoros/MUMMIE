#!/usr/bin/perl
use strict;

my %seeds;
open(IN,"all-seeds.fasta") || die;
while(<IN>) {
  chomp;
  if(/>/) {next}
  $seeds{$_}=1;
}
close(IN);

my @files=`ls predicted-seeds`;
foreach my $file (@files) {
  chomp $file;
  open(IN,"predicted-seeds/$file") || die;
  while(<IN>) {
    if(/>DNA/) {
      $_=<IN>;
      chomp;
      my $known=0+$seeds{$_};
      print "$known\t$_\n";
    }
  }
  close(IN);
}

