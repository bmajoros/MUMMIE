#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;
use TempFilename;

my $EXEC="/home/ohler/bmajoros/SCUMM";

my $name=ProgramName::get();
die "$name <*.hmm> <dir>\n" unless @ARGV==2;
my ($hmmFile,$dir)=@ARGV;

my $tmp1=TempFilename::generate();
my $tmp2=TempFilename::generate();

#open(OUT1,">tmp.predicted") || die;
#open(OUT2,">tmp.known") || die;
open(OUT1,">$tmp1") || die;
open(OUT2,">$tmp2") || die;
my @files=`ls $dir`;
foreach my $file (@files) {
  chomp $file;
  next unless $file=~/(.*).fastb$/;
  my $id=$1;
  my $cmd="$EXEC/parse -p $hmmFile $dir/$file";
  open(IN,"$cmd|") || die;
  while(<IN>) {print OUT1}
  close(IN);
  open(IN,"$dir/$id.path") || die "$dir/$id.path";
  while(<IN>) {print OUT2}
  close(IN);
}
close(OUT1);
close(OUT2);

system("$EXEC/eval.pl $tmp2 $tmp1");
system("rm $tmp1 $tmp2");


sub getProfile
  {
    my ($filename)=@_;
    open(IN,$filename) || die $filename;
    my $profile={};
    while(<IN>) {
      if(/(\d+)/){++$profile->{$1}}
    }
    close(IN);
    return $profile;
  }



