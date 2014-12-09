#!/usr/bin/perl
use strict;
use ProgramName;
use FastaReader;

my $name=ProgramName::get();
die "$name <chunks-dir> <seeds.fasta>\n" unless @ARGV==2;
my ($dir,$seedsFile)=@ARGV;

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
  next unless $file=~/([^\/]+)\.fastb/;
  my $substrate=$1;
  open(IN,"$dir/$file") || die $file;
  while(<IN>) {
    if(/>DNA/) {
      my $seed=<IN>;
      chomp $seed;
      my $ids=$IDs{$seed};
      if(!defined($ids)) {next}#{print STDERR "undefined: $seed ($seedsFile)\n";die $file}
      foreach my $id (@$ids) {
	print "$substrate\t$id\n";
	last;###
      }	
    }
  }
  close(IN);
}



