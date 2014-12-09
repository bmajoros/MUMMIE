


OBSOLETE



#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <indir-fastb> <indir-fold> <outdir> <dashUvalue>\n" unless @ARGV==4;
my ($fastbDir,$foldDir,$outdir,$U)=@ARGV;
my @fastbFiles=`ls $fastbDir`;
foreach my $fastb (@fastbFiles) {
  $fastb=~/([^\/]+)\.fastb/;
  my $stem=$1;
  my $fastbFile="$fastbDir/$stem.fastb";
  open(FASTB,$fastbFile) || die "can't open $fastbFile\n";
  my $foldFile="$foldDir/$stem.fold";
  open(FOLD,$foldFile) || die "can't open $foldFile\n";
  for(my $i=0 ; $i<202-$U+1 ; ++$i) { <FOLD> }
  my $outfile="$outdir/$stem.fastb";
  open(OUT,">$outfile") || die "can't write to file $outfile\n";
  while(<FASTB>) { print OUT }
  close(FASTB);
  print OUT "%unpaired\n";
  my @F;
  while(<FOLD>) {
    #if(/(\S+)\s+(\S+)/) { print OUT "$2\n" }
    if(/(\S+)\s+(\S+)/) { push @F,$2 }
  }
  my $FN=@F;
  for(my $i=0 ; $i<$FN ; ++$i) {
    my $end=$i+$U;
    if($end>$FN) {$end=$FN}
    my $sum=0;
    my $SS=0;
    for(my $j=$i ; $j<$end ; ++$j) {
      $sum+=$F[$j];
      ++$SS;
    }
    my $ave=$sum/$SS;
    print OUT "$ave\n";
  }
  close(OUT);
  close(FOLD);
}




