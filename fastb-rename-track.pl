#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;
use Fastb;

my $name=ProgramName::get();
die "$name <in.fastb> <old-name> <new-name>\n" unless @ARGV==3;
my ($filename,$oldName,$newName)=@ARGV;

my $fastb=new Fastb($filename);
my $track=$fastb->getTrackByName($oldName);
$track->rename($newName);
$fastb->save($filename);


