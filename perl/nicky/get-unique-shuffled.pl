#!/usr/bin/perl
use strict;
use FastaReader;
use ProgramName;

my $name=ProgramName::get();
die "$name <unshuffled.fasta> <shuffled.fasta>\n" unless @ARGV==2;
my ($unshuffledFile,$shuffledFile)=@ARGV;

my %keep;
my $reader=new FastaReader($unshuffledFile);
while(1) {
  my ($def,$seq)=$reader->nextSequence();
  last unless $def;
  $def=~/>(\S+)/ || die;
  my $id=$1;
  $keep{$id}=1;
}
$reader->close();
$reader=new FastaReader($shuffledFile);
while(1) {
  my ($def,$seq)=$reader->nextSequence();
  last unless $def;
  $def=~/>(\S+)_sh/ || die;
  my $id=$1;
  next unless $keep{$id};
  print "$def$seq\n";
}
$reader->close();



