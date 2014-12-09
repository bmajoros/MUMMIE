#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <graph-dir>\n" unless @ARGV==1;
my ($dir)=@ARGV;
my @files=`ls $dir`;

@files=("Mummie","targetscanCONS","targetscanALL","targetscanUNCONS");

foreach my $file (@files) {
  open(IN,"$dir/$file") || die $file;
  while(<IN>) {
    chomp;
    if(/(\S+)\s+(\S+)/) {
      my ($Sn,$Sp)=($1,$2);
      $Sn=~/(\S+):(\S+)/;
      $Sn=$1;
      $Sp=~/(\S+):(\S+)/;
      my ($mean,$sd)=($1,$2);

      ###
      #$sd/=sqrt(1000);
      ###

      my $top=$mean+$sd;
      my $bottom=$mean-$sd;
      print "$Sn\t$mean\n";
      print "$Sn\t$top\n";
      print "$Sn\t$bottom\n";
      print "$Sn\t$mean\n";
    }
    else {print}
  }
  close(IN);
  print "\n";
}


