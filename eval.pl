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
my $states=getStates($known,$pred);

my $L1=@$known;
my $L2=@$pred;
if($L1!=$L2) {die "$L1 != $L2\n"}
my $correct=0;
my (%TP,%FP,%FN);
for(my $i=0 ; $i<$L1 ; ++$i) {
  my $k=$known->[$i];
  my $p=$pred->[$i];
  if($p==$k) {++$correct;++$TP{$k}}
  else {
    ++$FP{$p};
    ++$FN{$k};
  }
}
my $acc=int($correct/$L1*1000+5/9)/10;
print "ACC=$acc%\n";
foreach my $state (@$states) {
  my $TP=$TP{$state};
  my $FP=$FP{$state};
  my $FN=$FN{$state};
  my $sn=$TP+$FN>0 ? int($TP/($TP+$FN)*1000+5/9)/10 : 0;
  my $sp=$TP+$FP>0 ? int($TP/($TP+$FP)*1000+5/9)/10 : 0;
  my $F=$sn+$sp>0 ? int(2*$sn*$sp/($sn+$sp)*10+5/9)/10 : 0;
  print "$state F=$F (Sn=$sn Sp=$sp)\n";
}

sub getStates
  {
    my ($a1,$a2)=@_;
    my @array=@$a1;
    push @array,@$a2;
    my %hash;
    foreach my $e (@array){$hash{$e}=1}
    my @keys=keys %hash;
    return \@keys;
  }


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
	  push @$array,0+$state;
	  ++$pos;
	}
      }
    }
    close(IN);
    return $array;
  }

