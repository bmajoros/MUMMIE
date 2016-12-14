#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <in.fastb> <out.bed> <trackname> <threshold> <min-len>\n"
  unless @ARGV==5;
my ($infile,$outfile,$wantTrack,$threshold,$MIN_LEN)=@ARGV;

$infile=~/([^\/]+)\.fastb/ || die;
my $substrate=$1;
open(OUT,">$outfile") || die;
open(IN,$infile) || die;
my $x=0;
my $echo=0;
my $inSite=0;
my ($begin,$end,$sum);
my $nextID=1;
while(<IN>) {
  if(/[%>]\s*(\S+)/) {
    if($1 eq $wantTrack) { $echo=1 } else { $echo=0 }
    $x=0;
  }
  elsif(/(\S+)/) {
    if($echo) {
      my $y=0+$_;
      if($y>=$threshold) {
	if(!$inSite) { $inSite=1; $begin=$x+1; $sum=0 }
	$sum+=$y
      }
      else {
	if($inSite) {
	  $inSite=0;
	  $end=$x+1;
	  if($end-$begin>=$MIN_LEN) {
	    my $score=$sum/($end-$begin);
	    my $id="elem$nextID";
	    ++$nextID;
	    print OUT "$substrate\t$begin\t$end\t$id\t$score\n";
	  }
	}
      }
    }
    ++$x;
  }
}
if($inSite) {
  $inSite=0;
  $end=$x+1;
  if($end-$begin>=$MIN_LEN) {
    my $score=$sum/($end-$begin);
    my $id="elem$nextID";
    ++$nextID;
    print OUT "$substrate\t$begin\t$end\t$id\t$score\n";
  }
}
close(IN);
close(OUT);



