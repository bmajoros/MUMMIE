#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;

my %lengths;
my @files=`ls fold`;
my $numFiles=@files;
for(my $i=0 ; $i<$numFiles ; ++$i) {
  my $file=$files[$i];
  chomp $file;
  my $out=`fasta-seq-length.pl fold/$file`;
  $out=~/(\d+)\s+bp/ || die $out;
  my $len=$1;
  $file=~/(.*)\.fasta/ || die $file;
  my $transcriptID=$1;
  $lengths{$transcriptID}=$len;
}

process("AGO-assembled");
process("PUM2-assembled");
process("QKI-assembled");


#=======================================================

sub process
  {
    my ($dir)=@_;
    my @files=`ls $dir`;
    my $numFiles=@files;
    for(my $i=0 ; $i<$numFiles ; ++$i) {
      my $file=$files[$i];
      my ($len,$transcriptID);
      open(IN,"$dir/$file") || die "$dir/$file";
      while(<IN>) {
	if(/length=(\d+)/) {$len=$1}
	elsif(/transcriptID=(\S+)/) {$transcriptID=$1}
      }	
      close(IN);
      my $L=$lengths{$transcriptID};
      if($L!=$len) {print "$transcriptID\t$len\t$L\n"}
    }
  }




