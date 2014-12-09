#!/usr/bin/perl
use strict;
use ProgramName;
use FastaReader;
use Fastb;
use GffReader;

my $name=ProgramName::get();
die "$name <motifs.fasta> <predictions.gff>\n" unless @ARGV==2;
my ($motifFile,$gffFile)=@ARGV;

my %motifs;
my $reader=new FastaReader($motifFile);
while(1) {
  my ($defline,$sequence)=$reader->nextSequence();
  last unless $defline;
  $defline=~/>(\S+)/ || die;
  $motifs{substr($sequence,0,8)}->{$1}=1;      # 8mer-m1
  $motifs{substr($sequence,0,7)."A"}->{$1}=1;  # 8mer-A1
  $motifs{substr($sequence,1,7)}->{$1}=1;      # 7mer-m1
  $motifs{substr($sequence,1,6)."A"}->{$1}=1;  # 7mer-A1
  $motifs{substr($sequence,0,7)}->{$1}=1;      # 7mer-m8
  $motifs{substr($sequence,0,6)}->{$1}=1;      # 6mer3-8
  $motifs{substr($sequence,1,6)}->{$1}=1;      # 6mer1-7
}

my %beats;
my $reader=new GffReader();
my $bySubstrate=$reader->hashBySubstrate($gffFile);
my @substrates=keys %$bySubstrate;
my $numSubstrates=@substrates;
for(my $i=0 ; $i<$numSubstrates ; ++$i) {
  my $substrate=$substrates[$i];
  my $features=$bySubstrate->{$substrate};
  my $groups=GffReader::groupByOverlaps($features);
  my $numGroups=@$groups;
  #print "$numGroups groups on $substrate\n";
  for(my $j=0 ; $j<$numGroups ; ++$j) {
    #print "group $j\n";
    my $group=$groups->[$j];
    my $groupSize=@$group;
    my %byMicroRNA;
    foreach my $feature (@$group) {
      my $seq;
      my $extra=$feature->{additionalFields};
      foreach my $field (@$extra) {
	#print "extra=\"$field\"\n";
	if($field=~/seq=(\S+);/) { $seq=$1 }
      }
      my $miRNAs=$motifs{$seq};
      #print "\"$seq\"\n";
      my @miRNAs=keys %$miRNAs;
      my $key=join(",",@miRNAs);
      push @{$byMicroRNA{$key}},$feature;
    }
    my @keys=keys %byMicroRNA;
    foreach my $key (@keys) {
      #print "$key:\n";
      my $features=$byMicroRNA{$key};
      my ($bestScore,$bestLength,%lengths);
      foreach my $feature (@$features) {
	#my $gff=$feature->toGff(); print "\t$gff";
	my $len=$feature->getLength();
	$lengths{$len}=1;
	my $score=$feature->getScore();
	if($score>$bestScore) {$bestScore=$score;$bestLength=$len}
      }
      my @lengths=keys %lengths;
      foreach my $len (@lengths) {
	next if $len==$bestLength;
	++$beats{$bestLength}->{$len};
      }
    }
  }
}
for(my $i=6 ; $i<=8 ; ++$i) {
  for(my $j=6 ; $j<=8 ; ++$j) {
    next if $i==$j;
    my $n=$beats{$i}->{$j};
    print "$i beats $j => $n times\n";
  }
}


