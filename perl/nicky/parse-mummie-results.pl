#!/usr/bin/perl
use strict;

my $echo;
open(IN,"Mummie.results.raw") || die;
while(<IN>) {
  if(/MUMMIE/) {$echo=1}
  #if($_=~/NAIVE/ || $_=~/eval-shuffle/) {$echo=0;print "\n"}
  if($_=~/NAIVE/) {$echo=0;print "\n"}
  next unless $echo;
  if(/(\S+):\s+(\S+)\s\+\/-\s+(\S+)/) {
    print "$2:$3\t";
  }
}
close(IN);

