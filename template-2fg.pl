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
die "$name <schema.txt> <order> <outfile>\n" unless @ARGV==3;
my ($schemaFile,$order,$outfile)=@ARGV;

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
4 states
schema:
$numDiscrete\t$numContinuous
";
foreach my $x (@discrete) {print OUT "$x\nACGT\n"}
foreach my $x (@continuous) {print OUT "$x\n"}
print OUT "
transitions:
4 4
0.0  0.001   0.000  0.00
1.0  0.990   0.100  0.10
0.0  0.0045  0.900  0.00
0.0  0.0045  0.000  0.90

state 1 emissions:
3       1    0    0
";
emission();
print OUT "
state 2 emissions:
3       0    1    0
";
emission();
print OUT "
state 3 emissions:
3       0    0    1
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
    if($numContinuous==1) {
print OUT "
1 0.00001
1        1
0.00001

1 1
1        1
0.1

1 1
1        1
0.1
";
}
    elsif($numContinuous==2) {
print OUT "
2 0.00001 0.5
2        2
0.00001  0.00001
0.00001  0.5

2 1      1
2        2
0.1      0.00001
0.00001  0.1

2 1      0.5
2        2
0.1      0.00001
0.00001  0.5
";
}
    else { die }
#    component($BG_MEAN,$BG_VAR,$BG_COV);
#    component($FG_MEAN,$FG_VAR,$FG_COV);
#    component($FG_MEAN,$FG_VAR,$FG_COV);
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
	#$chain.="$nmer\t-1.38629436111989\n";
	$chain.="$nmer\t$p\n";
      }
    }
    return $chain;
  }
