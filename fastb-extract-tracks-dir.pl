#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;

my $BIN=$ENV{"MUMMIE"};

my $name=ProgramName::get();
die "$name <indir> <outdir> <track1> <track2> ...\n" unless @ARGV>=3;
my $indir=shift @ARGV;
my $outdir=shift @ARGV;
my $tracks;
while(@ARGV) {
  my $track=shift @ARGV;
  $tracks.=" $track";
}

my @files=`ls $indir`;
foreach my $file (@files) {
  chomp $file;
  if($file=~/\.fastb$/) {
    system("$BIN/fastb-extract-tracks.pl $indir/$file $outdir/$file $tracks");
  }
}
