#!/usr/bin/perl
use strict;

#my $ids="pos1-ids.txt";
#my $ids="pos3-ids.txt";
my $ids="neg-ids.txt";

my %keep;
open(IN,$ids) || die $ids;
while(<IN>) {
  chomp;
  if(/(\S+)\s+(\S+)/) {
    $keep{$1}=1;
  }
}
close(IN);

open(IN,"targetscan.gff") || die;
while(<IN>) {
  chomp;
  my @fields=split/\t/,$_;
  if($keep{$fields[0]}) { print "$_\n" }
}
close(IN);



