#!/usr/bin/perl
use strict;
use GffReader;
use ProgramName;

my $MARGIN=15;
my $name=ProgramName::get();
die "$name <predictions.gff> <groups.gff>\n" unless @ARGV==2;
my ($predFile,$groupFile)=@ARGV;

for(my $iter=0 ; $iter<=100 ; ++$iter) {
  my $threshold=$iter/100;
  my $reader=new GffReader;
  my $predHash=$reader->hashBySubstrate($predFile);
  my $groupHash=$reader->hashBySubstrate($groupFile);
  my $predictions=[];
  my @substrates=keys %$predHash;
  for(my $i=0 ; $i<@substrates ; ++$i) {
    my $substrate=$substrates[$i];
    my $preds=$predHash->{$substrate};
    @$preds=sort {$a->getBegin() <=> $b->getBegin()} @$preds;
    my $n=@$preds;
    for(my $j=0 ; $j<$n ; ++$j) {
      my $this=$preds->[$j];
      if($this->getScore()<$threshold) {
	splice(@$preds,$j,1);
	--$n;
	--$j;
      }
    }
goto SKIP;
    for(my $j=0 ; $j+1<$n ; ++$j) {
      my $this=$preds->[$j];
      my $that=$preds->[$j+1];
      if($this->getBegin()==$that->getBegin() &&
	 $this->getEnd()==$that->getEnd()) {
	my $bestScore=$this->getScore()>$that->getScore() ? $this->getScore :
	  $that->getScore();
	$this->setScore($bestScore);
	splice(@$preds,$j+1,1);
	--$n;
	--$j;
      }
    }
  SKIP:
    #push @$predictions,@$preds;
    for(my $j=0 ; $j<$n ; ++$j) {
      my $prediction=$preds->[$j];
      my $add=join("",@{$prediction->{additionalFields}});
      my $keep=1;
      if($add=~/_sh_(\d+)/) {
	#if($1!=0) {$keep=0}
      }
      if($keep) {push @$predictions,$prediction }
    }
  }
  my $numPred=@$predictions;
  my ($TP,$FP);
  for(my $i=0 ; $i<$numPred ; ++$i) {
    my $prediction=$predictions->[$i];
    my $add=join("",@{$prediction->{additionalFields}});
    if($add=~/_sh/) {$prediction->{sh}=1}
    my $groups=$groupHash->{$prediction->getSubstrate()};
    my $n=$groups ? @$groups : 0;
    for(my $j=0 ; $j<$n ; ++$j) {
      my $group=$groups->[$j];
      if($prediction->overlapsOther($group)) {
	if($prediction->{sh}) {++$FP}
	else { $prediction->{mark}=$group->{mark}=1; ++$TP }
      }
    }
  }
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
  my $Sp=$hitGroups/$numPreds;
  my $F=2*$Sn*$Sp/($Sn+$Sp);
  print "threshold=$threshold Sn=$Sn Sp=$Sp F=$F hitGroups=$hitGroups numGroups=$numGroups numPreds=$numPreds\n";
}



