#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;

die "assemble-transcripts.pl <indir> <outdir>\n" unless @ARGV==2;
my ($indir,$outdir)=@ARGV;

my %transcripts;
my @files=`ls $indir`;
my $numFiles=@files;
for(my $i=0 ; $i<$numFiles ; ++$i) {
  my $file=$files[$i];
  $file=~/(\S+)_(\d+)_(\d+)\.fastb/ || die "can't parse filename: $file\n";
  my ($gene,$transcript,$exon)=($1,$2,$3);
  my $base="$gene\_$transcript";
  $transcripts{$base}++;
}
my @transcripts=keys %transcripts;
my $n=@transcripts;
for(my $i=0 ; $i<$n ;++$i) {
  my $base=$transcripts[$i];
  my $numExons=$transcripts{$base};
  my $outfile="$outdir/$base.fastb\n";
  my %tracks;
  my ($currentTrack,%trackTypes,%deflines);
  for(my $j=0 ; $j<$numExons ; ++$j) {
    my $infile="$indir/$base\_$j.fastb";
    open(IN,$infile) || die "can't open file $infile\n";
    while(<IN>) {
      if(/^([>%])(\S+)(.*)/) {
	my ($type,$id,$rest)=($1,$2,$3);
	$currentTrack=$id;
	$trackTypes{$id}=$type;
	$deflines{$id}=$_;
      }
      else {
	$tracks{$currentTrack}.=$_;
      }
    }
    close(IN);
  }
  open(OUT,">$outfile") || die "can't write to file $outfile\n";
  my @tracks=keys %tracks;
  foreach my $trackID (@tracks) {
    my $data=$tracks{$trackID};
    my $defline=$deflines{$trackID};
    if($defline=~/^%/) {
      $defline=~/^(.*\s+\/length=)\d+\s*$/ || die "can't parse: $defline";
      my $firstPart=$1;
      my $L=($data=~s/\n/\n/g);
      $defline="$firstPart$L\n";
    }
    else {
      $defline=~s/\s+\/begin=\d+\s+\/end=\d+//g || die "can't parse: $defline";
      chomp $defline;
      $defline.=" /numexons=$numExons\n";
    }
    print OUT "$defline$data";
  }
  close(OUT);
}


