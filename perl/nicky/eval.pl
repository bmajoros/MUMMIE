#!/usr/bin/perl
use strict;
use ProgramName;

# THIS PROGRAM EVALUATES ACCURACY FOR EXPERIMENT #1 (BHRF1/3)

my $name=ProgramName::get();
die "$name <pos1.txt> <pos2.txt <pos3.txt> <neg.txt> <Npos1> <Npos2> <Npos3> <Nneg>\n" unless @ARGV==8;
my ($file1,$file2,$file3,$file4,$N1,$N2,$N3,$N4)=@ARGV;

my ($junk1,$TP1)=analyze($file1,"ebv-mir-bhrf1-1");
my ($junk2,$TP2)=analyze($file2,"ebv-mir-bhrf1-2");
my ($junk3,$TP3)=analyze($file3,"ebv-mir-bhrf1-3");
my ($junk4,$FP)=analyze($file4,"ebv-mir-bhrf1-1","ebv-mir-bhrf1-2",
			"ebv-mir-bhrf1-3");

my $N=$N1+$N2+$N3+$N4;
my $TP=$TP1+$TP2+$TP3+$N3-$FP;
my $acc=$TP/$N;
my $TP=$TP1+$TP2+$TP3;
my $SNR=($TP/($N1+$N2+$N3)) / ($FP/$N4);
#my $SNR=$TP/$FP;
my $PPV=$TP/($TP+$FP);
my $FDR=$FP/($TP+$FP);
$acc=round($acc);
$SNR=round($SNR);
$PPV=round($PPV);
$FDR=round($FDR);
print "acc=$acc  SNR=$SNR  PPV=$PPV   FDR=$FDR\n";

sub round {
  my ($x)=@_;
  return int($x*1000+5/9)/1000;
}

sub analyze {
  my ($file,$mir,$mir2)=@_;
  my $TP=0;
  my $N=0;
  print "$file\n";
  my (%genes,%hits);
  open(IN,$file) || die "can't open file: $file\n";
  while(<IN>) {
    chomp;
    if(/(\S+)\s+(\S+)/) {
      my ($group,$predMir)=($1,$2);
      $group=~/^(\S+)-interval/ || die;
      my $gene=$1;
#print "gene=\"$gene\"\n";
      #my $gene=$group;
      $genes{$gene}=1;
      if($predMir eq $mir || $predMir eq $mir2) {$hits{$gene}=1}
    }
  }
  close(IN);
  my @keys=keys %genes;
  $N=@keys;
  @keys=keys %hits;
  $TP=@keys;
  print "$N $TP\n";
  return ($N,$TP);
}



