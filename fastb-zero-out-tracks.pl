#!/usr/bin/perl
use strict;
use ProgramName;
use Fastb;

my $name=ProgramName::get();
die "$name <in.fastb> <track,name,list> <out.fastb>\n" unless @ARGV==3;
my ($infile,$trackList,$outfile)=@ARGV;

my @tracks=split/,/,$trackList;
my $fastb=new Fastb($infile);
foreach my $trackName (@tracks) {
  my $track=$fastb->getTrackByName($trackName);
  my $L=$track->getLength();
  my $data=$track->getData();
  for(my $i=0 ; $i<$L ; ++$i) { $data->[$i]=0 }
}
$fastb->save($outfile);





