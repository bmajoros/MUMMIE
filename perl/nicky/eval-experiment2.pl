#!/usr/bin/perl
use strict;
use ProgramName;
use TempFilename;
use GffReader;

# THIS PROGRAM EVALUATES ACCURACY FOR EXPERIMENT #2 (shuffled/unshuffled)

my $name=ProgramName::get();
die "$name <predictions.gff>\n" unless @ARGV==1;
my ($gff)=@ARGV;
my $reader=new GffReader;
my $sites=$reader->loadGFF($gff);
my $N=@$sites;
my ($TP,$FP)=(0,0);
for(my $i=0 ; $i<$N ; ++$i) {
  my $site=$sites->[$i];
  my $fields=$site->{additionalFields};
  my $str=join(" ",@$fields);
  if($str=~/_sh/) { ++$FP }
  else { ++$TP }
  }
my $PPV=int($TP/($TP+$FP)*1000+5/9)/1000;
my $FDR=int($FP/($TP+$FP)*1000+5/9)/1000;
print "TP=$TP\tFP=$FP\tPPV=$PPV\tFDR=$FDR\n"




