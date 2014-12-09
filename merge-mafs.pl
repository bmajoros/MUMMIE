#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <indir> <outdir>\n" unless @ARGV==2;
my ($indir,$outdir)=@ARGV;

my @keep=("panTro2","gorGor1","ponAbe2","rheMac2","papHam1","calJac1","tarSyr1","micMur1","otoGar1");
my %keep;
foreach my $keep (@keep) {$keep{$keep}=1}

my @files=`ls $indir`;
@files=sort {$a <=> $b} @files;
my $n=@files;
my ($prev,@chunks);
for(my $i=0 ; $i<$n ; ++$i) {
  my $file=$files[$i];
  chomp $file;
  $file=~/(\S+)_(\d+)/ || die;
  my ($transcriptID,$chunkID)=($1,$2);
  if($transcriptID eq $prev) {
    push @chunks,$file;
  }
  else {
    @chunks=sort {$a=~/_(\d+)/;my $A=$1;$b=~/_(\d+)/;my $B=$1;$A<=>$B} @chunks;
    emit(\@chunks);
    @chunks=($file);
  }
  $prev=$transcriptID;
}
emit(\@chunks);


sub emit {
  my ($chunks)=@_;
  foreach my $chunk (@$chunks) {
    $chunk=~/(\S+)_(\d+)/ || die;
    my ($transcriptID,$chunkID)=($1,$2);
    my $outfile="$outdir/$transcriptID.fasta";
    my $file="$indir/$chunk";
    my %seqs;
    open(IN,$file) || die;
    while(<IN>) {
      my @fields=split/\s+/,$_;
      if($fields[0] eq "s") {
	my $species=$fields[1];
	my $seq=$fields[6];
	next unless($species=~/ENST/ || $keep{$species});
	$seqs{$species}.=$seq;
      }
    }
    close(IN);
    open(OUT,">$outfile") || die;
    my @keys=keys %seqs;
    foreach my $key (@keys) {
      my $seq=$seqs{$key};
      my $species=$key;
      if($species=~/ENST/) {$species="hg19"}
      print OUT ">$species\n$seq\n";
    }
    close(OUT);
  }
}






