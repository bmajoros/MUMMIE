#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <predictions.gff.txt> <scan.gff.txt>\n" unless @ARGV==2;
my ($file1,$file2)=@ARGV;

my $set1=load($file1);
my $set2=load($file2);
my @keys1=keys %$set1;
my @keys2=keys %$set2;
my $hits;
my $n1=@keys1;
my $n2=@keys2;
foreach my $key (@keys1) {
  if($set2->{$key}) {++$hits}
}
print "N1=$n1 N2=$n2 intersection=$hits\n";



sub load {
  my ($file)=@_;
  my $hash={};
  open(IN,$file) || die;
  while(<IN>) {
    if(/(\S+)-interval\d+\s+(\S+)/) {
      my ($substrate,$mir)=($1,$2);
      $hash->{$1}=$2;
    }
  }
  close(IN);
  return $hash;
}



