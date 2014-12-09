#!/usr/bin/perl
use strict;
use ProgramName;
use FastaReader;

my $name=ProgramName::get();
die "$name <chunks-dir> <seeds.fasta> <weights.txt>\n" unless @ARGV==3;
my ($dir,$seedsFile,$weightsFile)=@ARGV;

my %weights;
open(IN,$weightsFile) || die $weightsFile;
while(<IN>) {
  if(/(\S+)\s+(\S+)/) {
    my ($id,$weight)=($1,$2);
    $weights{$id}=$weight;
  }
}
close(IN);

my %IDs;
my $reader=new FastaReader($seedsFile);
while(1) {
  my ($defline,$sequence)=$reader->nextSequence();
  last unless $defline;
  $defline=~/>(\S+)/ || die "can't parse defline\n";
  my $id=$1;
  push @{$IDs{$sequence}},$id;
  push @{$IDs{substr($sequence,0,6)}},$id;
  push @{$IDs{substr($sequence,1,6)}},$id;
  push @{$IDs{substr($sequence,0,7)}},$id;
  push @{$IDs{substr($sequence,1,7)}},$id;
  push @{$IDs{substr($sequence,0,8)}},$id;
  push @{$IDs{substr($sequence,0,7)."A"}},$id;
  push @{$IDs{substr($sequence,1,6)."A"}},$id;
}
$reader->close();


my @files=`ls $dir`;
my $numFiles=@files;
for(my $i=0 ; $i<$numFiles ; ++$i) {
  my $file=$files[$i];
  chomp $file;
  next unless $file=~/\.fastb/;
  open(IN,"$dir/$file") || die $file;
  while(<IN>) {
    if(/>DNA/) {
      my $seed=<IN>;
      chomp $seed;
      my $ids=$IDs{$seed};
      if(!defined($ids)) {print STDERR "undefined: $seed\n";die $file}
      foreach my $id (@$ids) {
	my $weight=$weights{$id};
	print "$weight\n";
      }	
    }
  }
  close(IN);
}



