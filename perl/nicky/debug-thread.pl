#!/usr/bin/perl
use strict;
use GffReader;
use SummaryStats;
use ProgramName;
use TempFilename;
$|=1;

die "eval-targetscan-thread.pl <id> <seeds-file> <consY/N>\n" unless @ARGV==3;
my ($id,$seedsFile,$consYN)=@ARGV;
$id=~/(\S+)_(\S+)/ || die $id;
my ($sh,$threshold)=($1,$2);

my $CONSTRAIN_TO_GROUPS=1;
my $reader=new GffReader();
my $temp=TempFilename::generate();
print "tempfile=$temp\n";
my $groups=$reader->hashBySubstrate("groups.gff");

my $seeds=loadSeeds($seedsFile);
my ($shuffled,$unshuffled);
#-----------------------------------
if($consYN eq "Y") {
  extractShuffles("targetscan/sanity/shuffled-conserved.gff",$temp,$seeds);
  $shuffled=load($temp);
  $unshuffled=load("targetscan/sanity/unshuffled-conserved.gff");
}
else {
  extractShuffles("targetscan/sanity/shuffled.gff",$temp,$seeds);
  $shuffled=load($temp);
  $unshuffled=load("targetscan/sanity/unshuffled.gff");
}

my $filteredTrue=filterTP($unshuffled,$threshold);
my $filteredFalse=filterFP($shuffled,$threshold);
my $TP=@$filteredTrue;
my $FP=@$filteredFalse;
exit unless $FP;
my $filename=TempFilename::generate();
print "second tempfile=$filename\n";
open(OUT,">$filename") || die;
foreach my $p (@$filteredTrue) {
  print OUT $p->toGff();
print "keeping " . $p->toGff();
}
foreach my $p (@$filteredFalse) {
  print OUT $p->toGff();
print "keeping " . $p->toGff();
}
close(OUT);
open(IN,"eval-SnSp.pl $filename groups.gff |");
#open(IN,"eval-SnSp.pl $filename WT.clusters.gff |");
my (@SNR,@PPV,@FDR,@Sn,@Sp,@F);
while(<IN>) {
  if(/Sn=(\S+)\s+Sp=(\S+)\s+F=(\S+)/) {
    my ($sn,$sp,$f)=($1,$2,$3);
    push @Sn,$sn; push @Sp,$sp; push @F,$f;
  }
  elsif(/SNR=(\S+)\s+PPV=(\S+)\s+FDR=(\S+)/) {
    my ($snr,$ppv,$fdr)=($1,$2,$3);
    push @SNR,$snr; push @PPV,$ppv; push @FDR,$fdr;
  }
}
close(IN);

exit unless @SNR>0;
my ($mean,$stddev,$min,$max)=SummaryStats::roundedSummaryStats(\@SNR);
exit unless $mean>0;
print "$mean\t";
my ($mean,$stddev,$min,$max)=SummaryStats::roundedSummaryStats(\@PPV);
print "$mean\t";
my ($mean,$stddev,$min,$max)=SummaryStats::roundedSummaryStats(\@FDR);
print "$mean\t";

my ($mean,$stddev,$min,$max)=SummaryStats::roundedSummaryStats(\@Sn);
print "$mean\t";
my ($mean,$stddev,$min,$max)=SummaryStats::roundedSummaryStats(\@Sp);
print "$mean\t";
my ($mean,$stddev,$min,$max)=SummaryStats::roundedSummaryStats(\@F);
print "$mean\n";


sub filterTP {
  my ($sites,$threshold)=@_;
  my $selected=[];
  foreach my $site (@$sites) {
    if($site->getScore()<=$threshold) {
      my $found=0;
      my $G=$groups->{$site->getSubstrate()};
      #if(!$G) {print "no groups on " . $site->getSubstrate() . "\n"}
      if(!$G) {next}
      foreach my $group (@$G) {
	if($group->overlapsOther($site)) { $found=1; last }
      }
      if($found || !$CONSTRAIN_TO_GROUPS)
	{ push @$selected,$site }
    }
  }
  return $selected;
}
sub filterFP {
  my ($sites,$threshold)=@_;
  my $selected=[];
  foreach my $site (@$sites) {
    if($site->getScore()<=$threshold) {
      my $found=0;
      my $G=$groups->{$site->getSubstrate()};
      #if(!$G) {print "no groups on " . $site->getSubstrate() . "\n"}
      if(!$G) {next}
      foreach my $group (@$G) {
	if($group->overlapsOther($site)) { $found=1; last }
      }
      if($found || !$CONSTRAIN_TO_GROUPS)
	{ push @$selected,$site }
    }
  }
  return $selected;
}

sub load {
  my ($filename)=@_;
  my $features=$reader->loadGFF($filename);
  return $features;
}

sub loadSeeds {
  my ($filename)=@_;
  my $hash={};
  #my $count=0;
  open(IN,$filename) || die $filename;
  while(<IN>) {
    chomp;
    if(/>(\S+_sh_\S*)/) {
      $hash->{$1}=1;
      #++$count
    }
  }
  close(IN);
  #print "$count shuffled seeds\n";
  return $hash;
}

sub extractShuffles {
  my ($from,$to,$keep)=@_;
  open(IN,$from) || die $from;
  open(OUT,">$to") || die $to;
  while(<IN>) {
    if(/miRNA=([^;]+);/) {
      if($keep->{$1}) {
	print OUT;
	print "loading $_\n";
      }
    }
  }
  close(OUT);
  close(IN);
}




