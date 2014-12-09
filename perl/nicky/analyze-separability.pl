#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <pos.txt> <neg.txt>\n" unless @ARGV==2;
my ($file1,$file2)=@ARGV;

my (@A,@N);
load($file1,1);
load($file2,0);
@A=sort {$b->[0] <=> $a->[0]} @A;   # when higher scores are better
#@A=sort {$a->[0] <=> $b->[0]} @A;   # when lower scores are better
my $N=@A;
my $TP=0;
my $FP=0;
my $TN=$N[0];
my $FN=$N[1];
my $mod=int($N/100);
my ($AUC,$prevSn,$prevFPR);
for(my $i=0 ; $i<$N ; ++$i) {
  if($A[$i]->[1]==1) {
    ++$TP;
    --$FN;
  }
  else {
    ++$FP;
    --$TN;
  }
  if($i%$mod==0) {
    my $Sn=$TP/($TP+$FN);
    #my $Sp=$TP/($TP+$FP);
    #print "$Sn $Sp\n";
    my $FPR=$FP/($FP+$TN);
    print "$FPR $Sn\n";
    $AUC+=($FPR-$prevFPR)*($Sn+$prevSn)/2;
    $prevFPR=$FPR;
    $prevSn=$Sn;
  }
}
print STDERR "$AUC\n";

sub load {
  my ($filename,$label)=@_;
  open(IN,$filename) || die $filename;
  while(<IN>) {
    chomp;
    if(/(\S+)/) {
      push @A,[$1,$label];
      ++$N[$label];
    }
  }
  close(IN);
}
