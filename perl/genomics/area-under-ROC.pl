#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <points-file>\n" unless @ARGV==1;
my ($filename)=@ARGV;

my @points;
push @points,[1,1];
open(IN,$filename) || die "can't open file $filename\n";
while(<IN>) {
  if(/^\s*([\d\.e\+\-]+)\s+([\d\.e\+\-]+)/) {
    push @points,[$1,$2];
  }
}
close(IN);
push @points,[0,0];

my $area=0;
while(@points>1) {
  my $thisPoint=pop @points;
  my $nextPoint=$points[@points-1];
  my ($x1,$y1)=@$thisPoint;
  my ($x2,$y2)=@$nextPoint;
  #print "($x1,$y1) -> ($x2,$y2)  =  ";
  if($x1!=$x2 && $y1==$y2) {
    my $inc=$y1*($x2-$x1);
    $area+=$inc;
    #print "$inc";
  }
  #print "\n";
  #print "\t\t\t====> $area\n";
}
print "$area\n";


