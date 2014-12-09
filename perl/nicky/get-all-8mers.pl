#!/usr/bin/perl
use strict;
use ProgramName;
use FastaReader;

my $name=ProgramName::get();
die "$name <infile.fasta>\n" unless @ARGV==1;
my ($infile)=@ARGV;

my $reader=new FastaReader($infile);
while(1) {
  my ($def,$seq)=$reader->nextSequence();
  last unless $def;
  my $L=length($seq);
  for(my $i=0 ; $i<$L-8 ; ++$i) {
    my $sub=substr($seq,$i,8);
    print "$sub\n";
  }
}

