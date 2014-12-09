#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <model-filestem> <fastb-dir> <foreground-states>\n" unless @ARGV==3;
my ($MODEL,$FASTB_DIR,$FOREGROUND_STATES)=@ARGV;

my $OUTFILE="$MODEL.gff";
my $HMM_FILE="$MODEL.hmm";

system("rm $OUTFILE") if -e $OUTFILE;
my @files=`ls $FASTB_DIR`;
foreach my $file (@files) {
  chomp $file;
  next unless $file=~/\.fastb$/;
  system("MUMMIE/parse $HMM_FILE $FASTB_DIR/$file -p -g $FOREGROUND_STATES >> $OUTFILE");
}

