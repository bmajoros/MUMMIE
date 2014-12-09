#!/usr/bin/perl
use strict;
use GffReader;

my $reader=new GffReader;
#my $dashI=$reader->loadGFF("dashG.gff");
#my $hash=$reader->hashBySubstrate("dashI.gff");
my $dashI=$reader->loadGFF(" bulge1-groups.gff");
my $hash=$reader->hashBySubstrate("bulge1-groups-PD.gff");
my $n=@$dashI;
for(my $i=0 ; $i<$n ; ++$i) {
  my $feature=$dashI->[$i];
  my $substrate=$hash->{$feature->getSubstrate()};
  foreach my $other (@$substrate) {
    if($feature->getBegin()==$other->getBegin() &&
       $feature->getEnd()==$other->getEnd()) {
      my $s1=$feature->getScore();
      my $s2=$other->getScore();
      print "$s1\t$s2\n";
    }
  }
}



