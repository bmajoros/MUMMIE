#!/usr/bin/perl
use strict;
use GffReader;

#my $file1="pos1.txt";
#my $file2="pos3.txt";
#my $file3="neg.txt";

my $file1="pos1-proximal.txt";
my $file2="pos3-proximal.txt";
my $file3="neg-proximal.txt";

my $pos1Groups=loadGroups("pos1-groups.gff");
my $pos3Groups=loadGroups("pos3-groups.gff");
my $negGroups=loadGroups("neg-groups.gff");

my ($N1,$TP1)=analyze($file1,"ebv-mir-bhrf1-1",$pos1Groups);
my ($N2,$TP2)=analyze($file2,"ebv-mir-bhrf1-3",$pos3Groups);
my ($N3,$FP3)=analyze($file3,"ebv-mir-bhrf1-1",$negGroups);
my ($N3,$FP4)=analyze($file3,"ebv-mir-bhrf1-3",$negGroups);

my $N=$N1+$N2+$N3;
my $TP=$TP1+$TP2+$N3-$FP4;
my $acc=$TP/$N;
print "acc=$acc\n";

sub loadGroups {
  my ($filename)=@_;
  my $reader=new GffReader;
  my $features=$reader->loadGFF($filename);
  return $features;
}


sub analyze {
  my ($file,$mir)=@_;
  my $TP=0;
  my $N=0;
  my $numGroups=0;
  my (%genes,%hits);
  open(IN,$file) || die "can't open file: $file\n";
  while(<IN>) {
    chomp;
    if(/(\S+)\s+(\S+)/) {
      my ($group,$predMir)=($1,$2);
      $group=~/(\S+)-interval/ || die;
      my $gene=$1;
      $genes{$gene}=1;
      if($predMir eq $mir) {$hits{$gene}=1}
    }
  }
  close(IN);
  my @keys=keys %genes;
  $N=@keys;
  @keys=keys %hits;
  $TP=@keys;
  return ($N,$TP);
}



