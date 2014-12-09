#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <in.gff> <min-signal> <min-cons> <min-len>\n" unless @ARGV==4;
my ($infile,$minSig,$minCons,$minLen)=@ARGV;

open(IN,$infile) || die "can't open file: $infile\n";
LINE:
while(<IN>) {
  chomp;
  next if/^#/;
  my @fields=split/\s+/,$_;
  my $extra=$fields[8];
  my $len=$fields[4]-$fields[3]+1;
  next unless $len>=$minLen;
  my @extra=split/;/,$extra;
  foreach my $pair (@extra) {
    $pair=~/(\S+)=(\S+)/ || die "can't parse extra field: $pair\n";
    my $field=$1;
    my $value=$2;
    if($field=~/Signal/) { next LINE unless $value>=$minSig }
    elsif($field=~/phastcons/) { next LINE unless $value>=$minCons }
  }
  print "$_\n";
}
close(IN);

