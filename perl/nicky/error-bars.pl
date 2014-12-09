#!/usr/bin/perl
use strict;
use ProgramName;

while(<STDIN>) {
  chomp;
  if(/(\S+)\s+(\S+)/) {
    my ($Sn,$Sp)=($1,$2);
    $Sn=~/(\S+):(\S+)/;
    $Sn=$1;
    $Sp=~/(\S+):(\S+)/;
    my ($mean,$sd)=($1,$2);
    my $top=$mean+$sd;
    my $bottom=$mean-$sd;
    print "$Sn\t$mean\n";
    print "$Sn\t$top\n";
    print "$Sn\t$bottom\n";
    print "$Sn\t$mean\n";
  }
  else {print "\n"}
}
print "\n";



