#!/usr/bin/perl
use strict;
use GffReader;
use ProgramName;

my $name=ProgramName::get();
die "$name <predictions.gff> <seeds.fasta>\n" unless @ARGV==2;
my ($predictionFile,$seedsFile,$outfile)=@ARGV;

my %keep;
open(IN,$seedsFile) || die;
while(<IN>) {
  if(/>(\S+)/) {$keep{$1}=1}
}
close(IN);

open(IN,$predictionFile) || die;
while(<IN>) {
  chomp;
  next if(/#/);
  my @fields=split/\t/,$_;
  my $extra=$fields[8];
  $extra=~/miRNA=([^;]+);/;
  if($keep{$1}) {print "$_\n"}
}
close(IN);





