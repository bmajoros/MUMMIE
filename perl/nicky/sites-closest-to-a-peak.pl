#!/usr/bin/perl
use strict;
use ProgramName;
use GffReader;

my $name=ProgramName::get();
die "$name <peaks.gff> <scan.gff> <clusters.gff>\n" unless @ARGV==3 || @ARGV==2;
my ($peaksGFF,$scanGFF,$clustersGFF)=@ARGV;
my $IGNORE_CLUSTERS=($clustersGFF eq "");

my $reader=new GffReader;

my %clusters;
if($clustersGFF) {
  my $features=$reader->loadGFF($clustersGFF);
  my $n=@$features;
  for(my $i=0 ; $i<$n ; ++$i) {
    my $feature=$features->[$i];
    my $substrate=$feature->getSubstrate();
    my $begin=$feature->getBegin();
    my $end=$feature->getEnd();
    push @{$clusters{$substrate}},[$begin,$end];
  }
}

my (%peaks,%clusterByPeak);
my $features=$reader->loadGFF($peaksGFF);
my $n=@$features;
for(my $i=0 ; $i<$n ; ++$i) {
  my $feature=$features->[$i];
  my $substrate=$feature->getSubstrate();
  my $begin=$feature->getBegin();
  my $clusters=$clusters{$substrate};
  my $found=0;
  foreach my $cluster (@$clusters) {
    my ($from,$to)=@$cluster;
    if($begin>=$from && $begin<=$to) { $found=$cluster; last }
  }
  next unless $found || $clustersGFF eq "";
  push @{$peaks{$substrate}},$begin;
  $clusterByPeak{$begin}=$found;
}

my (@sites,$prevSubstrate);
my $features=$reader->loadGFF($scanGFF);
my $n=@$features;
for(my $i=0 ; $i<$n ; ++$i) {
  my $feature=$features->[$i];
  my $substrate=$feature->getSubstrate();
  my $begin=$feature->getBegin();
  if(defined($prevSubstrate) && $substrate ne $prevSubstrate) {
    pickSites($prevSubstrate);
    @sites=();
  }
  push @sites,$feature;
  $prevSubstrate=$substrate;
}
pickSites($prevSubstrate);


#=============================================================
sub pickSites
  {
    my ($substrate)=@_;
    my $peaks=$peaks{$substrate};
    if(!defined($peaks)) {return} #{ die "no peaks for $substrate" }
    my $numPeaks=@$peaks;
    my %sitesUsed;
    foreach my $peak (@$peaks) {
      my $cluster=$clusterByPeak{$peak};
      my ($bestSite,$bestDist,$bestScore);
      foreach my $site (@sites) {
	next if $sitesUsed{$site};
	if(!$IGNORE_CLUSTERS && isInCluster($site,$cluster)) {
	  print $site->toGff();
	  $sitesUsed{$site}=1;
	  next;
	}
	my $begin=$site->getBegin();
	my $end=$site->getEnd();
	my $score=$site->getScore();
	my $dist;
	if($peak<$begin) { $dist=$begin-$peak }
	elsif($peak>$end) { $dist=$peak-$end }
	else { $dist=0 }
	my $thisIn=isInCluster($site,$cluster);
	my $thatIn=$bestSite ? isInCluster($bestSite,$cluster) : 0;
	my $better=0;
	if(!defined($bestSite)) { $better=1 }
	elsif($thisIn && !$thatIn) { $better=1 }
	elsif(!$thisIn && $thatIn) { next }
	elsif($dist<$bestDist ||
	      $dist==$bestDist && $score>$bestScore) { $better=1 }
	if($better) {
	  $bestSite=$site;
	  $bestDist=$dist;
	  $bestScore=$score;
	}
      }
      next unless $IGNORE_CLUSTERS;
      next unless $bestSite;
      print $bestSite->toGff();
      $sitesUsed{$bestSite}=1;
    }
  }


sub isInCluster {
  my ($site,$cluster)=@_;
  if(!$cluster) { return 1 }
  my $siteBegin=$site->getBegin();
  my $siteEnd=$site->getEnd();
  my ($clusterBegin,$clusterEnd)=@$cluster;
  return $siteBegin<$clusterEnd && $clusterBegin<$siteEnd;
}
