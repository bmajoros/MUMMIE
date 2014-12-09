#!/usr/bin/perl
use strict;
use FastaReader;
use ProgramName;

my $name=ProgramName::get();
die "$name <in.fasta>\n" unless @ARGV==1;
my ($infile)=@ARGV;

my %seen;
my $reader=new FastaReader($infile);
while(1) {
  my ($def,$seq)=$reader->nextSequence();
  last unless $def;
  next if $seen{$seq};
  $seen{$seq}=1;
  print "$def$seq\n";
}

