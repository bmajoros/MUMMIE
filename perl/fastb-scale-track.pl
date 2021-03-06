#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use ProgramName;
use Fastb;

my $name=ProgramName::get();
die "$name <in-out.fastb> <track-name> <multiplier>\n" unless @ARGV==3;
my ($fastbFile,$trackName,$factor)=@ARGV;

my $fastb=new Fastb($fastbFile);
my $track=$fastb->getTrackByName($trackName);
my $data=$track->getData();
my $L=$track->getLength();
if(!$track->isContinuous()) { die "$trackName is non-numeric\n" }
for(my $i=0 ; $i<$L ; ++$i) { $data->[$i]*=$factor }
$fastb->save($fastbFile);



