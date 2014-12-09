#!/usr/bin/perl
use strict;
use Fastb;
use ProgramName;

my $name=ProgramName::get();
die "$name <fastb-dir> <min-readcount> <min-signal>\n" unless @ARGV==3;
my ($dir,$minReadCount,$minSignal)=@ARGV;

my ($goodGroups,$badGroups);
my @files=`ls $dir`;
my $numFiles=@files;
for(my $i=0 ; $i<$numFiles ; ++$i) {
  my $filename=$files[$i];
  chomp $filename;
  my $fastb=new Fastb("$dir/$filename");
  my $trReads=$fastb->getTrackByName("AGO-ReadCount");
  my $trSignal=$fastb->getTrackByName("AGO-ConversionPercent");
  my $trSignal=$fastb->getTrackByName("AGO-Signal");
  my $groups=$trReads->getNonzeroRegions();
  foreach my $group (@$groups) {
    my ($begin,$end)=@$group;
    my $readData=$trReads->getData();
    my $sigData=$trSignal->getData();
    my ($maxReads,$maxSignal);
    for(my $pos=$begin ; $pos<$end ; ++$pos) {
      my $x=$readData->[$pos];
      if($x>$maxReads) { $maxReads=$x }
      $x=$sigData->[$pos];
      if($x>$maxSignal) { $maxSignal=$x }
    }
    if($maxReads>$minReadCount) {
      if($maxSignal>$minSignal) {
	++$goodGroups;
      }
      else {
	++$badGroups;
      }
    }
    else {
    }
  }
}

my $total=$goodGroups+$badGroups;
my $percent=int(100*$goodGroups/$total+5/9);
print "$goodGroups groups with signal>$minSignal ($percent\%)\n";



