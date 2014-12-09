#!/usr/bin/perl
use strict;
use ProgramName;
use FastaReader;

my $name=ProgramName::get();
die "$name <scan-hits.gff> <fastb-dir> <track-name>\n" unless @ARGV==3;
my ($gffFile,$dir,$trackName)=@ARGV;

my %hits;
open(IN,$gffFile) || die;
while(<IN>) {
  chomp;
  my @fields=split/\s+/,$_;
  my ($gene,$id,$s,$begin,$end,$score)=@fields;
  push @{$hits{$gene}},[$begin-1,$end,$score,$id]
}
close(IN);

my %hist;
my @keys=keys %hits;
my $n=@keys;
my ($totalGroups,$totalLength,$numGroupsWithMultiplePeaks,$totalSites,$ebv,
    $groupLenSum,$fivePrime,$threePrime);
for(my $i=0 ; $i<$n ; ++$i) {
  my $key=$keys[$i];
  my $file="$dir/$key.fastb";
  open(IN,$file) || die "can't open $file";
  my (@array,$echo,$x);
  while(<IN>) {
    chomp;
    if(/[%>](\S+)/) {if($1 eq $trackName){$echo=1}}
    elsif($echo) {$array[$x++]=$_}
  }
  close(IN);
  my $hits=$hits{$key};
  my $numHits=@$hits;

  # find groups
  my $L=$x;
  $totalLength+=$L;
  my ($begin);
  for(my $pos=0 ; $pos<$L ; ++$pos) {
    my $y=$array[$pos];
    if($y>0 && ($pos==0 || $array[$pos-1]==0)) { $begin=$pos }
    elsif($y==0 && $pos>0 && $array[$pos-1]>0) {
      my $groupLen=$pos-$begin;
      $groupLenSum+=$groupLen;
      ++$totalGroups;

      # find peaks
      my $peaksThisGroup;
      for(my $x=$begin ; $x<$pos ; ++$x) {
	if(isPeak(\@array,$x,$L)) {
	  ++$peaksThisGroup;
	  for(my $j=0 ; $j<$numHits ; ++$j) {
	    #next if $hits->[$j]->[4];
	    my ($begin,$end,$score,$id)=@{$hits->[$j]};
	    my $siteLen=$end-$begin;
	    my $weight=$siteLen/21; # 21 = 6 + 7 +8

	    $weight=1; ###

	    next if($end<$x-20 || $begin>$x+20);
	    my $mid=($begin+$end)/2;
	    if($mid<$x) { $fivePrime+=$weight }
	    else { $threePrime+=$weight }
	    ++$totalSites;
	    if($id eq "ebv-mir-bhrf1-1" || $id eq "ebv-mir-bhrf1-3")
	      { ++$ebv }
	    #$hits->[$j]->[4]=1;
	  }
	}
      }
      if($peaksThisGroup>1) { ++$numGroupsWithMultiplePeaks }
    }
  }
}
my $pMultPeak=$numGroupsWithMultiplePeaks/$totalGroups;
print "P(mult peak)=$pMultPeak $numGroupsWithMultiplePeaks / $totalGroups\n";
my $probGroup=$totalGroups/($totalLength-$groupLenSum);
print "prob group = $probGroup $totalGroups / $totalLength\n";
my $pctFive=$fivePrime/($fivePrime+$threePrime);
my $pctThree=$threePrime/($fivePrime+$threePrime);
print "fiveprime = $pctFive $fivePrime, threeprime=$pctThree $threePrime $totalSites\n";
my $Pebv=$ebv/$totalSites;
print "P(ebv)=$Pebv\n";
my $aveGroupLen=$groupLenSum/$totalGroups;
print "ave group len = $aveGroupLen\n";

sub isPeak
  {
    my ($array,$pos,$L)=@_;
    my $y=$array->[$pos];
    if($pos+1<$L && $array->[$pos+1]>=$y) {return 0}
    if($pos>0 && $array->[$pos-1]>=$y) {return 0}
    return 1;
  }





