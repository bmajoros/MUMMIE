#!/usr/bin/perl
use strict;
use ProgramName;
$|=1;

my $name=ProgramName::get();
die "$name <model-filestem> <fastb-dir> <foreground-states> <outfile.gff>\n" unless @ARGV==4;
my ($MODEL,$FASTB_DIR,$FOREGROUND_STATES,$OUTFILE)=@ARGV;

my $NUM_JOBS=300;
#my $OUTFILE="$MODEL.gff";
my $HMM_FILE="$MODEL.hmm";

system("rm $OUTFILE") if -e $OUTFILE;
my @files=`ls $FASTB_DIR`;
my @cmds;
foreach my $file (@files) {
  chomp $file;
  next unless $file=~/\.fastb$/;
  my $cmd="MUMMIE/parse $HMM_FILE $FASTB_DIR/$file -p -g $FOREGROUND_STATES";
  push @cmds,$cmd;
}

system("rm parse-parallel.mpirun");
system("rm q/*");
system("rm parallel/*");
my $numCmds=@cmds;
my $cmdsPerJob=int($numCmds/$NUM_JOBS);
if($cmdsPerJob<1) {$cmdsPerJob=1}
my $nextCmdIndex=0;
for(my $i=0 ; $i<$NUM_JOBS ; ++$i) {
  my $qFile="q/$i.q";
  open(OUT,">$qFile") || die;
  print OUT 
"#!/bin/tcsh
#
#\$ -S /bin/tcsh -cwd
#\$ -o /home/ohler/bmajoros/MUMMIE/paper/peak-calling/AGO/nicky/parse-parallel.mpirun -j y
#\$ -l highprio
#\$ -N parse$i
cd /home/ohler/bmajoros/MUMMIE/paper/peak-calling/AGO/nicky
";
  for(my $j=0 ; $j<$cmdsPerJob && $nextCmdIndex<$numCmds ; ++$j) {
    my $cmd=$cmds[$nextCmdIndex++];
    print OUT "$cmd > parallel/$nextCmdIndex.gff\n";
  }
  close(OUT);
  last unless $nextCmdIndex<$numCmds;
  system("qsub $qFile >& /dev/null");
}
while(1) {
  my @qs=`qstat -u bmajoros`;
  my $n=@qs;
  last unless $n>0;
  print "waiting for $n jobs...\n";
  sleep(30);
}
system("cat parallel/*.gff > $OUTFILE");



