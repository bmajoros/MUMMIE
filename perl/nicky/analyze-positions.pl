#!/usr/bin/perl
use strict;
use GffReader;
use SummaryStats;

my $peaksFile="peaks.gff";
my $predFile="tmp.1";
my $reader=new GffReader();
my $predictions=$reader->loadGFF($predFile);
my $peaksHash=$reader->hashBySubstrate($peaksFile);

my %byMiRNA;
my $numPred=@$predictions;
for(my $i=0 ; $i<$numPred ; ++$i) {
  my $pred=$predictions->[$i];
  my $peaks=$peaksHash->{$pred->getSubstrate()};
  my ($bestPeak,$bestDist);
  foreach my $peak (@$peaks) {
    my $dist=abs($peak->getBegin()-$pred->getBegin());
    if(!defined($bestPeak) || $dist<$bestDist) {
      $bestPeak=$peak;
      $bestDist=$dist;
    }
  }
  next unless $bestPeak;
  my $pos=$pred->getBegin()-$bestPeak->getBegin();
  print "$pos\n";
  my %miRNAs;
  my $add=$pred->{additionalFields};
  foreach my $field (@$add) {
    if($field=~/miRNA=([^;]+)/) {$miRNAs{$1}=1}
  }
  my @miRNAs=keys %miRNAs;
  foreach my $miRNA (@miRNAs) {push @{$byMiRNA{$miRNA}},$pos}
}
my @keys=keys %byMiRNA;
foreach my $miRNA (@keys) {
  my $positions=$byMiRNA{$miRNA};
  my ($mean,$stddev,$min,$max)=SummaryStats::roundedSummaryStats($positions);
  $mean=int($mean+($mean>0?5/9:-5/9));
  $stddev=int($stddev+($stddev>0?5/9:-5/9));
  $min=int($min+($min>0?5/9:-5/9));
  $max=int($max+($max>0?5/9:-5/9));
  my $N=@$positions;
  #print "$miRNA:\t$mean +/- $stddev\t$min\-$max   N=$N\n";
}




