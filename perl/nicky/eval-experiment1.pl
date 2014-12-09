#!/usr/bin/perl
use strict;
use ProgramName;
use TempFilename;

# THIS PROGRAM EVALUATES ACCURACY FOR EXPERIMENT #1 (BHRF1/3)

my $G="G"; #"G";
#my $SCHEMA="WT.schema";
my $SCHEMA="WT-cons.schema";
my $SEEDS="nicky-top50-seeds.fasta";

my $name=ProgramName::get();
die "$name <pos1.gff> <pos2.gff> <pos3.gff> <neg.gff>\n" unless @ARGV==4;
my ($file1,$file2,$file3, $file4)=@ARGV;

my $tempDir=TempFilename::generate();
system("mkdir $tempDir");

process($file1,"pos1$G","$file1.txt");
process($file2,"pos2$G","$file2.txt");
process($file3,"pos3$G","$file3.txt");
process($file4,"neg$G","$file4.txt");
my $N1=0+`ls pos1$G | wc -l`;
my $N2=0+`ls pos2$G | wc -l`;
my $N3=0+`ls pos3$G | wc -l`;
my $N4=0+`ls neg$G | wc -l`;
system("eval.pl $file1.txt $file2.txt $file3.txt $file4.txt $N1 $N2 $N3 $N4");

system("rm -r $tempDir");

sub process {
  my ($gff,$chunks,$outFile)=@_;
  system("rm $tempDir/*");
  print "extracting $gff\n";
  system("MUMMIE/extract-gff-intervals $gff $chunks $SCHEMA $tempDir");
  print "done extracting\n";
  system("identify-seeds.pl $tempDir $SEEDS > $outFile");
}




