#!/usr/bin/perl
use strict;

  my $filename="seeds/0.seeds";
  my $hash={};
  open(IN,$filename) || die $filename;
  while(<IN>) {
    chomp;
    #if(/>(\S+_sh_\S*)/) {
    if(/>(\S+)/) {
      $hash->{$1}=1;
    }
  }
  close(IN);
my @keys=keys %$hash;
my $n=@keys;
print "$n\n";


filterSeeds("15percent/targetscan/shuffled.gff","tmp.tmp",$hash);

sub filterSeeds {
  my ($from,$to,$keep)=@_;
  open(IN,$from) || die $from;
  open(OUT,">$to") || die $to;
  while(<IN>) {
    if(/miRNA=([^;]+);/) {
      if($keep->{$1}) {print OUT}
    }
  }
  close(OUT);
  close(IN);
}



