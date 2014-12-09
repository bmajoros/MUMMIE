#!/usr/bin/perl
use strict;
use GffReader;
use ProgramName;

my $MARGIN=15;
my $name=ProgramName::get();
die "$name <predictions.gff> <scan.gff> <groups.gff>\n" unless @ARGV==3;
my ($predFile,$scanFile,$groupFile)=@ARGV;

my $reader=new GffReader;
my $scanHash=$reader->hashBySubstrate($scanFile);

for(my $iter=0 ; $iter<=100 ; ++$iter) {
  my $threshold=$iter/100;
  my $reader=new GffReader;
  my $predHash=$reader->hashBySubstrate($predFile);
  my $groupHash=$reader->hashBySubstrate($groupFile);
  my @substrates=keys %$predHash;

  foreach my $substrate (@substrates) {
    my $groups=$groupHash->{$substrate};
    my $scans=$scanHash->{$substrate};
    foreach my $scan (@$scans) {
      next unless $scan->getScore()>7; ##### USING 8MERS ONLY #####
      foreach my $group (@$groups) {
	if($scan->overlapsOther($group)) {
	  $group->{scanHit}=1;
	}
      }
    }
  }

  my $predictions=[];
  for(my $i=0 ; $i<@substrates ; ++$i) {
    my $substrate=$substrates[$i];
    my $preds=$predHash->{$substrate};
    #@$preds=sort {$a->getBegin() <=> $b->getBegin()} @$preds;
    @$preds=sort {
      $a->getBegin()<$b->getBegin() ? -1 :
	($a->getBegin()>$b->getBegin() ? 1 :
	 ($a->getEnd()<$b->getEnd() ? -1 :
          ($a->getEnd()>$b->getEnd() ? 1 : 0)
	 )
	)
      } @$preds;
    my $n=@$preds;
    for(my $j=0 ; $j<$n ; ++$j) {
      my $this=$preds->[$j];
      if($this->getScore()<$threshold) {
	splice(@$preds,$j,1);
	--$n;
	--$j;
      }
    }
    foreach my $prediction (@$preds) {
      my $add=join("",@{$prediction->{additionalFields}});
      if($add=~/_sh_/) {$prediction->{sh}=1}
    }
#goto SKIP;
    for(my $j=0 ; $j+1<$n ; ++$j) {
      my $first=$preds->[$j];
      my $last=$j;
      for(my $k=$j+1 ; $k<$n ; ++$k) {
	my $next=$preds->[$k];
	if($first->getBegin()!=$next->getBegin() ||
	   $first->getEnd()!=$next->getEnd()) { last }
	$last=$k;
      }
      next unless $last>$j;
      my ($numShuffled,$numNonShuffled);
      for(my $k=$j ; $k<=$last ; ++$k) 
	{ if($preds->[$k]->{sh}) {++$numShuffled} else {++$numNonShuffled} }
      if($numShuffled>0 && $numNonShuffled>0) {
	my $numDel=$last-$j+1;
	splice(@$preds,$j,$numDel);
	$n-=$numDel;
	--$j;
      }
      else {
	my $numDel=$last-$j;
	splice(@$preds,$j+1,$numDel);
	$n-=$numDel;
	--$j;
      }
    }
#  SKIP:
    for(my $j=0 ; $j<$n ; ++$j) {
      my $prediction=$preds->[$j];
      push @$predictions,$prediction;
    }
  }
  my $numPred=@$predictions;
  my ($TP,$FP);
  for(my $i=0 ; $i<$numPred ; ++$i) {
    my $prediction=$predictions->[$i];
    my $groups=$groupHash->{$prediction->getSubstrate()};
    my $n=$groups ? @$groups : 0;
    for(my $j=0 ; $j<$n ; ++$j) {
      my $group=$groups->[$j];
      if($prediction->overlapsOther($group) && !$group->{scanHit}) {
	$prediction->{mark}=1;
	if($prediction->{sh}) {++$FP}
	else { $group->{mark}=1; ++$TP }
      }
    }
  }
  if($TP && !$FP) {die "$TP $FP"}
  next unless $FP>0;
  my $SNR=$TP/$FP;
  my $PPV=$TP/($TP+$FP);
  my $FDR=$FP/($TP+$FP);
  print "threshold=$threshold SNR=$SNR PPV=$PPV FDR=$FDR TP=$TP FP=$FP\n";
  my @allGroups;
  my @substrates=keys %$groupHash;
  my $numSubstrates=@substrates;
  for(my $i=0 ;$i<$numSubstrates ; ++$i) {
    my $groups=$groupHash->{$substrates[$i]};
    push @allGroups,@$groups;
  }
  my $hitGroups=0;
  my $numGroups=@allGroups;
  my $hitPreds=0;
  my $numPreds=@$predictions;
  for(my $i=0 ; $i<$numGroups ; ++$i) {
    my $group=$allGroups[$i];
    if($group->{mark}) { ++$hitGroups }
  }
  for(my $i=0 ; $i<$numPreds ; ++$i) {
    my $pred=$predictions->[$i];
    if($pred->{mark}) { ++$hitPreds }
  }
  
  my $Sn=$hitGroups/$numGroups;
  my $Sp=$hitGroups/$hitPreds;
  my $F=2*$Sn*$Sp/($Sn+$Sp);
  print "threshold=$threshold Sn=$Sn Sp=$Sp F=$F hitGroups=$hitGroups numGroups=$numGroups numPreds=$numPreds\n";
}



