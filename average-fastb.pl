#!/usr/bin/perl
use strict;
use Fastb;
use ProgramName;

my $name=ProgramName::get();
die "$name <dir> <out.fastb>\n" unless @ARGV==2;
my ($dir,$outfile)=@ARGV;

my @files=`ls $dir`;
my $numFiles=@files;
my ($master,$L);
for(my $i=0 ; $i<$numFiles ; ++$i) {
  my $file=$files[$i];
  my $fastb=new Fastb("$dir/$file");
  if(!$master) { $master=$fastb; $L=$master->getLength(); next }
  my $numTracks=$fastb->numTracks();
  for(my $j=0 ; $j<$numTracks ; ++$j) {
    my $track=$fastb->getIthTrack($j);
    next unless $track->getType() eq "continuous";
    my $name=$track->getID();
    my $masterTrack=$master->getTrackByName($name);
    my $masterData=$masterTrack->getData();
    my $data=$track->getData();
    for(my $pos=0 ; $pos<$L ; ++$pos) {
      my $x=$data->[$pos];
      $masterData->[$pos]+=$data->[$pos];
    }
  }
  undef $fastb;
}
my $numTracks=$master->numTracks();
for(my $j=0 ; $j<$numTracks ; ++$j) {
  my $track=$master->getIthTrack($j);
  next unless $track->getType() eq "continuous";
  my $data=$track->getData();
  for(my $pos=0 ; $pos<$L ; ++$pos) {
    $data->[$pos]/=$numFiles;
  }
}
$master->save($outfile);



