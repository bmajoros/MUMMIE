#!/usr/bin/perl
use strict;

my %seeds;
open(IN,"/home/ohler/bmajoros/nicky/15percent/mirbase/top100.fasta") || die;
while(<IN>) {
  chomp;
  if(/>(\S+)/) {
    $_=<IN>;
    chomp;
    $seeds{$1}=$_;
  }
}

my %seen;
while(<STDIN>) {
  chomp;
  my @fields=split/\t/,$_;
  my ($substrate,$binding,$site,$begin,$end)=@fields;
  my $key="$substrate $begin $end";
  next if $seen{$key};
  #if(/(.*)seq=(\S+);type=8mer-A1;miRNA=(\S+);/) {
  if(/(.*)seq=(\S+);type=8mer-m1;miRNA=(\S+);/) {
    my ($stuff,$seq,$mir)=($1,$2,$3);
    my $seed=$seeds{$mir};
    my $expected=substr($seed,7,1);
    my $observed=substr($seq,7,1);
    #die unless $observed eq "A";
    if($observed eq "A" && $expected eq "A") {
      #print "${stuff}seq=$seq;type=8mer-m1;miRNA=$mir;\n";
      print "${stuff}seq=$seq;type=8mer-A1;miRNA=$mir;\n";
      $seen{$key}=1;
      next;
    }
  }
  print "$_\n";
  $seen{$key}=1;
}


