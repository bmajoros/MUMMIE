#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <main-track> <fastb-dir>\n" unless @ARGV==2;
my ($mainTrackName,$fastbDir)=@ARGV;

my @files=`ls $fastbDir`;
my $numFiles=@files;
my (%sums,%antisums,$N,$antiN);
for(my $i=0 ; $i<$numFiles ; ++$i) {
  my $file=$files[$i];
  open(IN,"$fastbDir/$file") || die "can't open file: $fastbDir/$file\n";
  my $echo=0;
  my ($currentTrack,%tracks);
  while(<IN>) {
    chomp;
    if(/>/) {$echo=0}
    elsif(/%(\S+)/) {$echo=1; $currentTrack=$1}
    elsif($echo) {push(@{$tracks{$currentTrack}},0+$_)}
  }
  close(IN);
  my @trackNames=keys %tracks;
  my $numTracks=@trackNames;
  my $mainTrack=$tracks{$mainTrackName};
  die "no track named $mainTrackName\n" unless $mainTrack;
  my $L=@$mainTrack;
  for(my $pos=0 ; $pos<$L ; ++$pos) {
    my $main=$mainTrack->[$pos];
    for(my $t=0 ; $t<$numTracks ; ++$t) {
      my $trackName=$trackNames[$t];
      next if $trackName eq $mainTrackName;
      my $track=$tracks{$trackName};
      my $val=$track->[$pos];
      if($main>0) { $sums{$trackName}+=$val; ++$N }
      else { $antisums{$trackName}+=$val; ++$antiN }
    }
  }
}
my @trackNames=keys %sums;
foreach my $trackName (@trackNames) {
  my $ave=$sums{$trackName}/$N;
  my $ave2=$antisums{$trackName}/$antiN;
  print "$trackName\t$ave\t$ave2\n";
}



