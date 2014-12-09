#!/usr/bin/perl
use strict;
use POSIX;

my $usage="$0 <infile.ps> <outfile.ps>";
die "$usage\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

my (@pages,$header,$currentPage);
open(IN,$infile) || die "Can't open file: \"$infile\"\n";
while(<IN>)
  {
    if(/^\%\%Page:\s+(\d+)\s+(\d+)\s*$/)
      {
	if($1==1) {$header=$currentPage}
	else {push @pages,$currentPage}
	$currentPage="";
      }
    else {$currentPage.=$_;}
  }
if(length($currentPage)>0) {push @pages,$currentPage}
close(IN);

open(OUT,">$outfile") || die "Can't creat file: \"$outfile\"\n";
print OUT $header;
my $n=@pages;
my $sheets=POSIX::ceil($n/4);
for(my $i=0 ; $i<$sheets ; ++$i)
  {
    my $first=$i*4;
    printPage($pages[$first+3]);
    printPage($pages[$first]);
    printPage($pages[$first+1]);
    printPage($pages[$first+2]);
  }
print OUT "showpage\n";
close(OUT);

sub printPage
  {
    my ($page)=@_;
    if(defined($page)) {print OUT $page}
    else {print OUT "showpage\n"}
  }
