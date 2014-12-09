#!/usr/bin/perl
use strict;
use ProgramName;

my $FASTB_DIR="fastb-WT-cons-longest";#"fastb-WT";

my $POSDIR1="pos1";
my $POSDIR3="pos3";
my $NEGDIR="neg";
system("mkdir $POSDIR1") unless -e $POSDIR1;
system("mkdir $POSDIR3") unless -e $POSDIR3;
system("mkdir $NEGDIR") unless -e $NEGDIR;

my ($genes1,$UTRs1)=loadScan("mir1-scan.gff");
my ($genes3,$UTRs3)=loadScan("mir3-scan.gff");
my $wtList=getGenes($FASTB_DIR);
my $diffD1=getGenes("WT-D1");
my $diffD3=getGenes("WT-D3");
my $sub1=subtract($wtList,$diffD1);
my $sub13=subtract($sub1,$diffD3);
my @negUTRs;
foreach my $gene (@$sub13) {
  #push @negUTRs,"$gene\_0";
  push @negUTRs,$gene;
}

copy($UTRs1,$POSDIR1);
copy($UTRs3,$POSDIR3);
copy(\@negUTRs,$NEGDIR);

my $n1=@$sub1;
my $n13=@$sub13;
my $nWT=@$wtList;
my $ndiff1=@$diffD1;
my $ndiff2=@$diffD3;

#print "sub1=$n1 sub13=$n13 wt=$nWT diff1=$ndiff1 diff2=$ndiff2\n";

#emit($sub13);

sub copy {
  my ($IDs,$dir)=@_;
  foreach my $id (@$IDs) {
    my $cmd="cp $FASTB_DIR/$id.fastb $dir";
    #print "$cmd\n";
    system($cmd);
  }
}

sub emit {
  my ($list)=@_;
  my $n=@$list;
  for(my $i=0 ; $i<$n ; ++$i) {
    my $elem=$list->[$i];
    print "$elem\_0\n";
  }
}

sub subtract {
  my ($A,$B)=@_;
  my $nA=@$A;
  my $nB=@$B;
  my $R=[];
  for(my $i=0 ; $i<$nA ; ++$i) {
    my $a=$A->[$i];
    my $found=0;
    for(my $j=0 ; $j<$nB ; ++$j) {
      my $b=$B->[$j];
      if($a eq $b) { $found=1; last }
    }
    if(!$found) { push @$R,$a }
  }
  return $R;
}

sub getGenes {
  my ($dir)=@_;
  my @files=`ls $dir`;
  my $genes=[];
  my $n=@files;
  my %hash;
  for(my $i=0 ; $i<$n ; ++$i) {
    my $file=$files[$i];
    $file=~/([^\/]+)\.fastb/ || die $file;
    my $iso=$1;
    #$iso=~/(\S+)_\d+/ || die;
    my $gene=$1;
    $hash{$gene}=1;
  }
  @$genes=keys %hash;
  return $genes;
}

sub loadScan {
  my ($filename)=@_;
  my %byGene;
  open(IN,$filename) || die $filename;
  while(<IN>) {
    chomp;
    my @fields=split/\s+/,$_;
    my $field1=$fields[0];
    my $iso=$field1;
    if($field1=~/(\S+)-/) { $iso=$1 }
    $iso=~/(\S+)_\d+/ || die;
    #my $gene=$1;
    my $gene=$iso;
    $byGene{$gene}=$iso; # pick an isoform arbitrarily
  }
  close(IN);
  my $values=[]; @$values=values %byGene;
  my $keys=[]; @$keys=keys %byGene;
  return ($keys,$values);
}



