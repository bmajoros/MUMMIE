#!/usr/bin/perl

$argc=@ARGV;
die "minus.pl <set1-file> <set2-file>\n" unless $argc==2;
($filename1,$filename2)=@ARGV;

open(FILE,$filename2) || die "Can't open $filename2\n";
while(<FILE>)
  {
    if($_=~/^\s*(\S.*\S)\s*$/ || $_=~/^\s*(\S+)\s*$/)
      {
	$hash{$1}=1;
      }
  }
close(FILE);

open(FILE,$filename1) || die "Can't open $filename1\n";
while(<FILE>)
  {
    if($_=~/^\s*(\S.*\S)\s*$/ || $_=~/^\s*(\S+)\s*$/)
      {
	print "$1\n" unless $hash{$1};
      }
  }
close(FILE);



