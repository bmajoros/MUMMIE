#!/usr/bin/perl
use strict;
use ProgramName;
use GffReader;

my $name=ProgramName::get();
die "$name <predictions.gff> <groups-dir> <out-dir>\n" unless @ARGV==3;
my ($gffFile,$groupsDir,$outDir)=@ARGV;

my $reader=new GffReader();
my $hash=$reader->hashBySubstrate($gffFile);
#my @keys=keys %$hash;

my @files=`ls $groupsDir`;
my $numFiles=@files;
for(my $i=0 ; $i<$numFiles ; ++$i) {
  my $filename=$files[$i];
  chomp $filename;
  $filename=~/([^\/]+)\.fastb/ || die;
  my $gene=$1;
  next if $hash->{$gene};
  system("cp $groupsDir/$filename $outDir");
}


