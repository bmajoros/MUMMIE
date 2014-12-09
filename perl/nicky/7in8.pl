#!/usr/bin/perl
use strict;
use FastaReader;
use ProgramName;

my @alpha=('A','C','G','T');

my $name=ProgramName::get();
die "$name <in.fasta>\n" unless @ARGV==1;
my ($infile)=@ARGV;

my $reader=new FastaReader($infile);
while(1) {
  my ($def,$seq)=$reader->nextSequence();
  last unless $def;
  chomp $def;
  my $L=length($seq);
  for(my $pos=0 ; $pos<$L ; ++$pos) {
    my $ref=substr($seq,$pos,1);
    foreach my $mut (@alpha) {
      next if $mut eq $ref;
      my $new=$seq;
      substr($new,$pos,1)=$mut;
      print "$def\_$pos$mut\n$new\n";
    }
  }
}
