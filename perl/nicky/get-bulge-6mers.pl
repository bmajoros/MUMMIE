#!/usr/bin/perl
use strict;
use ProgramName;
use FastaReader;

my $name=ProgramName::get();
die "$name <8mer-seeds.fasta> <out1.fasta> <out2.fasta>\n" unless @ARGV==3;
my ($seedFile,$outfile1,$outfile2)=@ARGV;

open(OUT1,">$outfile1") || die "can't write to file: $outfile1\n";
open(OUT2,">$outfile2") || die "can't write to file: $outfile2\n";
my $reader=new FastaReader($seedFile);
while(1) {
  my ($defline,$sequence)=$reader->nextSequence();
  last unless $defline;
  my $six=substr($sequence,1,6);
  my $firstTwo=substr($six,0,2);
  my $lastFour=substr($six,2,4);
  my $firstThree=substr($six,0,3);
  my $lastThree=substr($six,3,3);
  my $second=substr($six,1,1);
  my $third=substr($six,2,1);
  my $form1="$firstTwo$second$lastFour";
  my $form2="$firstThree$third$lastThree";
  print OUT1 "$defline$form1\n";
  print OUT2 "$defline$form2\n";
}
close(OUT1);
close(OUT2);



