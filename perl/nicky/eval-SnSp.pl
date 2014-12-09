#!/usr/bin/perl
use strict;
use GffReader;
use ProgramName;

my $name=ProgramName::get();
die "$name <predictions.gff> <groups.gff>\n" unless @ARGV==2;
my ($predFile,$groupFile)=@ARGV;

my $USE_GROUPS=1;

my $reader=new GffReader;
my $predHash=$reader->hashBySubstrate($predFile);
my $groupHash=$reader->hashBySubstrate($groupFile);
my $predictions=[];
my @substrates=keys %$predHash;
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
#SKIP:
  for(my $j=0 ; $j<$n ; ++$j) {
    my $prediction=$preds->[$j];
    push @$predictions,$prediction;
  }
}

my $numPred=@$predictions;
my $FP_outsideGroups=0; my $TP_outsideGroups=0;
my $TP=0; my $FP=0;
for(my $i=0 ; $i<$numPred ; ++$i) {
  my $prediction=$predictions->[$i];
  my $groups=$groupHash->{$prediction->getSubstrate()};
  my $n=$groups ? @$groups : 0;
  for(my $j=0 ; $j<$n ; ++$j) {
    my $group=$groups->[$j];
    if($prediction->overlapsOther($group)) {
      push @{$group->{sites}},$prediction;
    }
  }
}
for(my $i=0 ; $i<$numPred ; ++$i) {
  my $prediction=$predictions->[$i];
  my $groups=$groupHash->{$prediction->getSubstrate()};
  my $n=$groups ? @$groups : 0;
  for(my $j=0 ; $j<$n ; ++$j) {
    my $group=$groups->[$j];
    if($prediction->overlapsOther($group)) {
      $prediction->{mark}=1;
      if($prediction->{sh}) {++$FP}
      else { $group->{mark}=1; ++$TP }
    }
  }
  if(!$prediction->{mark}) {
    if($prediction->{sh}) { ++$FP_outsideGroups }
    else { ++$TP_outsideGroups }
  }
}

my $SNR;
if(!$USE_GROUPS) { $TP+=$TP_outsideGroups; $FP+=$FP_outsideGroups }
if($FP==0) { exit }
elsif($TP==0) { $SNR=0 }
else { $SNR=$TP/$FP }
my $PPV=$TP+$FP>0 ? $TP/($TP+$FP) : 0;
my $FDR=$TP+$FP>0 ? $FP/($TP+$FP) : 0;
print "SNR=$SNR PPV=$PPV FDR=$FDR TP=$TP FP=$FP\n";
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
my $Sp=$hitGroups/($USE_GROUPS ? $hitPreds : $numPreds);
my $F=$Sn+$Sp>0 ? 2*$Sn*$Sp/($Sn+$Sp) : 0;
print "Sn=$Sn Sp=$Sp F=$F hitGroups=$hitGroups numGroups=$numGroups numPreds=$numPreds\n";

