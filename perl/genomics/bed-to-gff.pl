#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <in.bed> <source-ID> <out.gff>\n" unless @ARGV==3;
my ($infile,$source,$outfile)=@ARGV;

open(OUT,">$outfile") || die "Can't write to file: $outfile\n";
open(IN,$infile) || die "Can't open file: $infile\n";
while(<IN>) {
  my @fields=split/\s+/,$_;
  next unless @fields>=6;
  my ($substrate,$begin,$end,$id,$zero,$strand)=@fields;
  my $extra;
  if(@fields>6) {
    for(my $i=6 ; $i<@fields ; ++$i) {
      $extra.="\t";
      $extra.=$fields[$i];
    }
  }
  print OUT "$substrate\t$source\t$id\t$begin\t$end\t$strand\t$zero$extra\n";
}
close(IN);
close(OUT);


# BED: chr10   82085865        82085894        NM_138812       0       -
# GFF: chr3R region   1   4526526  4526790  +  0

