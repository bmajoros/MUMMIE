#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;

my $BASEDIR="/home/ohler/bmajoros/SCUMM";

my $name=ProgramName::get();
die "$name <in.hmm> <outdir> <numseqs> <minlen> <maxlen>\n" unless @ARGV==5;
my ($hmmFile,$outdir,$numSeqs,$minLen,$maxLen)=@ARGV;

if(-e $outdir) {die "$outdir already exists\n"}
system("mkdir -p $outdir");

for(my $i=0 ; $i<$numSeqs ; ++$i) {
  my $id=$i+1;
  my $seed=abs(int(rand(1000000000)));
  my $cmd="$BASEDIR/sample -r $seed -m $minLen -M $maxLen $hmmFile $outdir/$id.fastb $outdir/$id.path";
  system($cmd);
}



