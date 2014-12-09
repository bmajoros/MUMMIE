#!/usr/bin/perl
use strict;
use GffReader;
use SummaryStats;
use ProgramName;
$|=1;




OBSOLETE





exit





my $name=ProgramName::get();
die "$name <seeds-dir-in> <seeds-dir-out> <conserved-Y/N>\n" unless @ARGV==3;
my ($inSeeds,$outSeeds,$consYN)=@ARGV;

system("rm -r $outSeeds") if -e $outSeeds;
system("mkdir $outSeeds");
for(my $sh=0 ; $sh<100 ; ++$sh) {
  system("cat $inSeeds/$sh.training-seeds.fasta | grep '>' > $outSeeds/$sh.seeds");
}

my $CONSTRAIN_TO_GROUPS=1;
my $reader=new GffReader();
#my $groups=$reader->hashBySubstrate("WT.clusters.gff");
my $groups=$reader->hashBySubstrate("groups.gff");
for(my $threshold=-1 ; $threshold<=0.4 ; $threshold+=0.1) {
  my (@SNR,@PPV,@FDR,@Sn,@Sp,@F);
  for(my $sh=0 ; $sh<100 ; ++$sh) {
    my $seeds=loadSeeds("$outSeeds/$sh.seeds");
    my ($shuffled,$unshuffled);
    #-----------------------------------
    if($consYN eq "Y") {
      extractShuffles("targetscan/group-conserved-top50.gff","tmp.1",$seeds);
      $shuffled=load("tmp.1");
      $unshuffled=load("targetscan/unshuffled-conserved-top50.gff");
    }
    else {
      #===
      extractShuffles("targetscan/group-top50.gff","tmp.1",$seeds);
      $shuffled=load("tmp.1");
      $unshuffled=load("targetscan/unshuffled-top50.gff");
      #system("cat targetscan/group-top50.gff | grep _sh_$sh\\; > tmp.1");
      #my $shuffled=load("tmp.1");
      #my $unshuffled=load("targetscan/unshuffled-top50.gff");
    }
    #-----------------------------------

    my $filteredTrue=filterTP($unshuffled,$threshold);
    my $filteredFalse=filterFP($shuffled,$threshold);
    my $TP=@$filteredTrue;
    my $FP=@$filteredFalse;
    next unless $FP;
    my $filename="tmp.tgsn";
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
  }
  next unless @SNR>0 && @Sn>0;
  my ($mean,$stddev,$min,$max)=SummaryStats::roundedSummaryStats(\@SNR);
  next unless $mean>0;
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
}

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
      if($keep->{$1}) {print OUT}
    }
  }
  close(OUT);
  close(IN);
}


