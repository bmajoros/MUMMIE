#!/usr/bin/perl
use strict;
use ProgramName;
use FastaReader;
use Translation;

my $name=ProgramName::get();
die "$name <in.fasta>\n" unless @ARGV==1;
my ($infile)=@ARGV;

my (@seqs,%seen);
my $reader=new FastaReader($infile);
while(1) {
  my ($def,$seq)=$reader->nextSequence();
  last unless $def;
  chomp $def; chomp $seq;
  push @seqs,[$def,$seq];
  $seen{$seq}=1;
}
my $n=@seqs;
my $hits=0;
my %emitted;
for(my $i=0 ; $i<$n ; ++$i) {
  my $pair=$seqs[$i];
  my ($def,$seq)=@$pair;
  #my $rev=reverse $seq;
  my $rev=Translation::reverseComplement(\$seq);
  if($seen{$rev}) { ++$hits; next }
  if($emitted{$rev}) {next}
  print "$def\_rev\n$rev\n";
  $emitted{$rev}=1;
}
print STDERR "$hits COLLISIONS\n";


