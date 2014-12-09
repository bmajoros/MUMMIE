#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <in.txt>\n" unless @ARGV==1;
my ($infile)=@ARGV;

my %ids;
open(IN,"transcript-ids.txt") || die;
while(<IN>) {
  chomp;
  my @fields=split/\s+/,$_;
  my ($gene,$transcript)=@fields;
  $ids{$gene}=$transcript;
}
close(IN);

open(IN,$infile) || die;
while(<IN>) {
  chomp;
  my $transcript=$ids{$_};
  print "$transcript\n";
}
close(IN);


