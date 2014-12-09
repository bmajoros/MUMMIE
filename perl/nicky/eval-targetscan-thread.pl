#!/usr/bin/perl
use strict;
use GffReader;
use SummaryStats;
use ProgramName;
use TempFilename;
$|=1;

#my $WANT_NON_CONS=1;

die "eval-targetscan-thread.pl <id> <seeds-file> <consY/N> <NONCONS:0/1>\n" unless @ARGV==4;
my ($id,$seedsFile,$consYN,$WANT_NON_CONS)=@ARGV;
$id=~/(\S+)_(\S+)/ || die $id;
my ($sh,$threshold)=($1,$2);

my $CONSTRAIN_TO_GROUPS=1; # 0=no, 1=yes

my $reader=new GffReader();
my $temp=TempFilename::generate();
#$temp="/nfs/igsp/ohler_fc/ohlerlab/bmajoros/COMMIE/temp/$temp";
my $groups=$reader->hashBySubstrate("groups.gff");

my $seeds=loadSeeds($seedsFile);
my ($shuffled,$unshuffled);
#-----------------------------------
my $BASEDIR="15percent/targetscan";
if($consYN eq "Y") {
  filterSeeds("$BASEDIR/shuffled-conserved.gff",
	      $temp,$seeds);
  $shuffled=load($temp);
  filterSeeds("$BASEDIR/unshuffled-conserved.gff",
	      $temp,$seeds);
  $unshuffled=load($temp);
}
else {
  filterSeeds($WANT_NON_CONS ? "$BASEDIR/shuffled-NONconserved.gff" :
	      "$BASEDIR/shuffled.gff",$temp,$seeds);
  $shuffled=load($temp);
  filterSeeds($WANT_NON_CONS ? "$BASEDIR/unshuffled-NONconserved.gff" :
	      "$BASEDIR/unshuffled.gff",$temp,$seeds);
  $unshuffled=load($temp);
}
system("rm $temp");

my $filteredTrue=filterGroupAndThreshold($unshuffled,$threshold);
my $filteredFalse=filterGroupAndThreshold($shuffled,$threshold);
my $TP=@$filteredTrue;
my $FP=@$filteredFalse;
exit unless $FP; ###
my $filename="targetscan/temp/$id.gff";#TempFilename::generate();
open(OUT,">$filename") || die;
foreach my $p (@$filteredTrue) {
  print OUT $p->toGff();
}
foreach my $p (@$filteredFalse) {
  print OUT $p->toGff();
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
#system("rm $filename");

exit unless @SNR>0;
my ($mean,$stddev,$min,$max)=SummaryStats::roundedSummaryStats(\@SNR);
#exit unless $mean>0;
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


sub filterGroupAndThreshold {
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
  open(IN,$filename) || die $filename;
  while(<IN>) {
    chomp;
    #if(/>(\S+_sh_\S*)/) {
    if(/>(\S+)/) {
      $hash->{$1}=1;
    }
  }
  close(IN);
  return $hash;
}

sub filterSeeds {
  my ($from,$to,$keep)=@_;
  open(IN,$from) || die $from;
  open(OUT,">$to") || die $to;
  while(<IN>) {
    if(/miRNA=([^;]+);/) {
      if($keep->{$1}) {print OUT}
    }
  }
  close(OUT);
  close(IN);
}




