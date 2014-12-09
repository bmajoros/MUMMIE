#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <sequence> <out.hmm>\n" unless @ARGV==2;
my ($seq,$outfile)=@ARGV;

my $L=length($seq);
my $numStates=$L+1;
open(OUT,">$outfile") || die "can't write to file: $outfile\n";
print OUT "$numStates states
schema:
1       0
DNA
ACGT

transitions:
$numStates $numStates
";
for(my $to=0 ; $to<$numStates ; ++$to) {
  for(my $from=0 ; $from<$numStates ; ++$from) {
    my $P=0;
    if($from+1==$to) { $P=1 }
    elsif($from==$numStates-1 && $to==0) { $P=1 }
    print OUT "$P\t";
  }
  print OUT "\n";
}
print OUT "\n";
for(my $q=1 ; $q<$numStates ; ++$q) {
  print OUT "state $q emissions:\n";
  print OUT "0 order\nalphabet:\nACGT\n";
  my $s=substr($seq,$q-1,1);
  foreach my $a ("A","C","G","T") {
    my $P="-inf";
    if($a eq $s) { $P=1 }
    print OUT "$a\t$P\n";
  }
}
close(OUT);


