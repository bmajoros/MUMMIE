#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <infile.gff>\n" unless @ARGV==1;
my ($infile)=@ARGV;

my $LEN=7;
my $MIN_CONS=0.8;

open(IN,$infile) || die;
while(<IN>) {
  chomp;
  my @fields=split/\s+/,$_;
  my $len=$fields[4]-$fields[3]+1;
  next unless $len==$LEN;
  $fields[8]=~/Signal=([^;]+);/ || die;
  my $ago=$1;
  next unless $ago>0;
  $fields[8]=~/phastcons=([^;]+);/ || die;
  my $cons=$1;
  next unless $cons>=$MIN_CONS;
  print "$_\n";
}
close(IN);

