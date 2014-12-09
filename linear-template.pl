#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;
use NgramIterator;

my $FG_MEAN=0.5;
my $BG_MEAN=0.0001;
my $BG_VAR=0.9;
my $BG_COV=0.01;
my $FG_VAR=0.2;
my $FG_COV=0.1;

my $name=ProgramName::get();
die "$name <schema.txt> <order> <motif-length> <outfile>\n" unless @ARGV==4;
my ($schemaFile,$order,$motifLen,$outfile)=@ARGV;
my $NUM_STATES=4+$motifLen;

my (@continuous,@discrete);
open(IN,$schemaFile) || die "Can't open $schemaFile\n";
while(<IN>) {
  if(/(\S+)\s*:\s*continuous/) {
    push @continuous,$1;
  }
  elsif(/(\S+)\s*:\s*discrete/) {
    push @discrete,$1;
  }
}
close(IN);
my $numDiscrete=@discrete;
my $numContinuous=@continuous;

open(OUT,">$outfile") || die "Can't write to file $outfile\n";
print OUT "
$NUM_STATES states
schema:
$numDiscrete\t$numContinuous
";
foreach my $x (@discrete) {print OUT "$x\nACGT\n"}
foreach my $x (@continuous) {print OUT "$x\n"}
print OUT "
transitions:
$NUM_STATES $NUM_STATES
";

my @T;
$T[0]->[1]=1;
$T[1]->[0]=0.001;
$T[1]->[1]=0.99;
$T[1]->[2]=0.009;
$T[2]->[1]=0.05;
$T[2]->[2]=0.9;
$T[2]->[3]=0.05;
for(my $q=3 ; $q<$NUM_STATES-1 ; ++$q) {
  $T[$q]->[$q+1]=1;
}
$T[$NUM_STATES-1]->[2]=0.99;
$T[$NUM_STATES-1]->[3]=0.01;
for(my $to=0 ; $to<$NUM_STATES ; ++$to) {
  for(my $from=0 ; $from<$NUM_STATES ; ++$from) {
    my $p=0+$T[$from]->[$to];
    print OUT "$p\t";
  }
  print OUT "\n";
}

print OUT "
state 1 emissions:
2 1 0
";
emission();
for(my $q=2 ; $q<$NUM_STATES ; ++$q) {
  print OUT "
state $q emissions:
2 0 1
";
  emission();
}
close(OUT);


#================================================================

sub component
  {
    my ($mean,$var,$cov)=@_;
    print OUT "\n$numContinuous";
    for(my $i=0 ; $i<$numContinuous ; ++$i) { print OUT "\t$mean" }
    print OUT "\n$numContinuous\t$numContinuous\n";
    for(my $i=0 ; $i<$numContinuous ; ++$i) {
      for(my $j=0 ; $j<$numContinuous ; ++$j) {
	if($i==$j) {print OUT "$var\t"}
	else {print OUT "$cov\t"}
      }
      print OUT "\n";
    }
  }

sub emission
  {
    component($BG_MEAN,$BG_VAR,$BG_COV);
    component($FG_MEAN,$FG_VAR,$FG_COV);
    my $CHAIN=loadChain($order);
    print OUT "\n$CHAIN";
  }

sub loadChain
  {
    my ($order)=@_;
    my $chain="$order order\nalphabet:\nACGT\n";
    for(my $n=1 ; $n<=$order+1 ; ++$n) {
      my $ngramIterator=new NgramIterator("ACGT",$n);
      my $nmer;
      while($nmer=$ngramIterator->nextString()) {
	my $p=log(rand(0.9998)+0.0001);
	$chain.="$nmer\t$p\n";
      }
    }
    return $chain;
  }
