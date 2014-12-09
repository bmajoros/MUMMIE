#!/usr/bin/perl
use strict;
use GffReader;
use Fastb;

my $chunksDir="fastb-WT-cons-longest";
my $peaksFile="WT-peaks.gff";
my $MARGIN=30;

my $reader=new GffReader();
my $peaksByChunk=$reader->hashBySubstrate($peaksFile);

my @counts;
my $numCounts=0;
my @chunkIDs=keys %$peaksByChunk;
my $numChunks=@chunkIDs;
for(my $i=0 ; $i<$numChunks ; ++$i) {
  my $chunkID=$chunkIDs[$i];
  my $peaks=$peaksByChunk->{$chunkID};
  my $fastbFilename="$chunksDir/$chunkID.fastb";
  next unless -e $fastbFilename;
  my $fastb=new Fastb($fastbFilename);
  my $track=$fastb->getTrackByName("phastcons");
  die $fastbFilename unless $track;
  my $L=$track->getLength();
  my $data=$track->getData();
  foreach my $record (@$peaks) {
    my $peak=$record->getBegin();
    my $left=$peak-$MARGIN;
    my $right=$peak+$MARGIN;
    next unless $left>=0 && $right<$L;
    for(my $pos=$left ; $pos<=$right ; ++$pos) {
      my $cons=$data->[$pos];
      $counts[$pos-$left]+=$cons;
      ++$numCounts;
    }
  }
}
my $L=@counts;
for(my $pos=0 ; $pos<$L ; ++$pos) {
  my $cons=$counts[$pos]/$numCounts;
  print "$pos\t$cons\n";
}
