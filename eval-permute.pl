#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <known> <predicted>\n" unless @ARGV==2;
my ($knownFile,$predFile)=@ARGV;
my $known=load($knownFile);
my $pred=load($predFile);

my $L1=@$known;
my $L2=@$pred;
if($L1!=$L2) {die "$L1 != $L2\n"}
my $correct=0;
for(my $i=0 ; $i<$L1 ; ++$i) {
  if($known->[$i]==$pred->[$i]) {++$correct}
}
my $acc=int($correct/$L1*1000+5/9)/10;
print "$acc%\n";

sub load
  {
    my ($filename)=@_;
    my $array=[];
    my $pos=0;
    open(IN,$filename) || die $filename;
    while(<IN>) {
      if(/>/) {}
      elsif(/(\S+)/) {
	my @points=split/\s+/,$_;
	foreach my $state (@points){
	  push @$array,$state;
	  ++$pos;
	}
      }
    }
    close(IN);
    return $array;
  }

