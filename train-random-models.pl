#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;

my $BASE="/home/ohler/bmajoros/SCUMM/test/rand";
my $BIN="/home/ohler/bmajoros/SCUMM";

for(my $i=0 ; $i<100 ; ++$i) {
  my $id=$i+1;
  my $modelName="rand$id";
  system("cd $BASE ; $BIN/random-HMM 4 1.0 8 random.schema $modelName.hmm ; sleep 1");
  writeQFile($modelName);
}



#=====================================================================
sub writeQFile
  {
    my ($modelName)=@_;
    my $qfile="$BASE/$modelName.q";
    my $initHmm="$BASE/$modelName.hmm";
    my $finalHmm="$BASE/$modelName-trained.hmm";
    open(OUT,">$qfile") || die $qfile;
    print OUT
"#!/bin/tcsh
#
#\$ -S /bin/tcsh -cwd
#\$ -o $BASE/mpirun.$modelName -j y
#\$ -pe threaded 8
#\$ -l highprio
#\$ -N EM$modelName
cd $BASE
$BIN/baum-welch $initHmm complete.tgf train 100 $finalHmm -c 8 -R >> $modelName.stdout
";
    close(OUT);
  }
#=====================================================================

