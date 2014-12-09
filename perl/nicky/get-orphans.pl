#!/usr/bin/perl
use strict;

my %hit;
open(IN,"hit-groups.txt") || die;
while(<IN>) {
  chomp;
  if(/(\S+)/) {$hit{$1}=1}
}
close(IN);

my @files=`ls groups2`;
foreach my $file (@files) {
  chomp $file;
  $file=~/(.*).fastb/ || next;
  my $stem=$1;
  if($hit{$stem}) {
    system("MUMMIE/fastb-to-fasta.pl groups2/$file non-orphans/$stem.fasta");
  }
  else {
    system("MUMMIE/fastb-to-fasta.pl groups2/$file orphans/$stem.fasta");
  }
}


