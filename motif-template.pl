#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;
use NgramIterator;

my $FG_MEAN=0.1;
my $BG_MEAN=0.0001;
my $BG_VAR=0.9;
my $BG_COV=0.01;
my $FG_VAR=0.2;
my $FG_COV=0.1;

my $name=ProgramName::get();
die "$name <schema.txt> <order> <outfile>\n" unless @ARGV==3;
my ($schemaFile,$order,$outfile)=@ARGV;
my $CHAIN=loadChain($order);

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
3 states
schema:
$numDiscrete\t$numContinuous
";
foreach my $x (@discrete) {print OUT "$x\nACGT\n"}
foreach my $x (@continuous) {print OUT "$x\n"}
print OUT "
transitions:
3 3
0       0.001   0
1       0.99    0.1
0       0.01    0.9

state 1 emissions:
2       0.99    0.01
";
emission();
print OUT "
state 2 emissions:
2       0.01    0.99
";
emission();
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
    print OUT "\n$CHAIN";
  }

sub loadChain2
  {
    my $chain;
    open(IN,"/home/ohler/bmajoros/COMMIE/random-chain.txt") || die;
    for(my $i=0 ; $i<3 ; ++$i) { $chain.=<IN> }
    while(<IN>) {
      #$chain.=$_;
      if(/(\S+)\s+(\S+)/) {$chain.="$1\t-1.38629436111989\n"}
    }
    close(IN);
    return $chain;
  }


sub loadChain
  {
    my ($order)=@_;
    my $chain="$order order\nalphabet:\nACGT\n";
    for(my $n=1 ; $n<=$order+1 ; ++$n) {
      my $ngramIterator=new NgramIterator("ACGT",$n);
      my $nmer;
      while($nmer=$ngramIterator->nextString()) {
	$chain.="$nmer\t-1.38629436111989\n"
      }
    }
    return $chain;
  }
