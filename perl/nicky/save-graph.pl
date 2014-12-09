#!/usr/bin/perl
use strict;

my $highest=0;
my @files=`ls graphs`;
foreach my $file (@files) {
  chomp $file;
  if($file=~/^(\d+)$/) {
    if($1>$highest) {$highest=$1}
  }
}
++$highest;
cmd("mkdir graphs/$highest");
cmd("cp Mummie targetscanCONS targetscanALL targetscanUNCONS graphs/$highest");
print "saved to graphs/$highest\n";

sub cmd {
  my ($cmd)=@_;
  #print "$cmd\n";
  system($cmd);
}



