#!/usr/bin/perl
use strict;
use FastaReader;

my ($numInexact,$numFound,$numTaxa);
my $bhrf1="TCAGGTTA";
open(IN,"cat 8mer-transcripts.txt 7mer-transcripts.txt |") || die;
while(<IN>) {
  chomp;
  if(/(\S+)/) {
    my $transcriptID=$1;
    my $filename="alignments-assembled/$transcriptID.fasta";
    if(-e $filename) {
      my ($foundExact,$foundInexact,$numTaxa)=(0,0,$numTaxa);
      my $reader=new FastaReader($filename);
      while(1) {
	my ($def,$seq)=$reader->nextSequence();
	last unless $def;
	++$numTaxa;
	$def=~/>(\S+)/ || die $def;
	my $taxon=$1;
	next if($taxon=~/^hg\d+/);
	$seq=~s/[-.]//g;
	if($seq=~/(TCAGGTTA)/) {++$foundExact;print "$1\n";}
	elsif(inexact($seq,$bhrf1)) {++$foundInexact}
      }
      $reader->close();
      #print "$foundExact exact\t$foundInexact inexact\t$numTaxa taxa\n";
    }
  }
}
close(IN);
#print "$numFound have exact matches\n";
#print "$numInexact have inexact matches\n";

sub inexact {
  my ($seq,$pattern)=@_;
  my $L=length($pattern);
  my $len=length($seq);
  for(my $i=0 ; $i<$len-$L ; ++$i) {
    my $nmer=substr($seq,$i,$L);
    my $numMismatch=countMismatches($nmer,$pattern);
    if($numMismatch==1) {
      print "$nmer\n";
      return 1
    }
  }
  return 0;
}

sub countMismatches {
  my ($seq,$pattern)=@_;
  my $L=length($seq);
  my $n=0;
  for(my $i=0 ; $i<$L ; ++$i) {
    if(substr($seq,$i,1) ne substr($pattern,$i,1)) { ++$n }
  }
  #print "$seq $pattern $n\n";
  return $n;
}




