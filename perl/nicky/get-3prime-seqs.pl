#!/usr/bin/perl
use strict;
use ProgramName;
use Fastb;
use FastaReader;
use GffReader;
use Translation;
use SmithWaterman;

my $name=ProgramName::get();
die "$name <fastb-dir> <sites.gff> <mature.fasta> <out-dir>\n" unless @ARGV==4;
my ($inDir,$gffFile,$matureFile,$outDir)=@ARGV;

my $aligner=new SmithWaterman(".","dna.matrix");
my $matureHash=FastaReader::readAll($matureFile); # hash : id->sequence
my $reader=new GffReader;
my $sitesHash=$reader->hashBySubstrate($gffFile);
my @files=`ls $inDir`;
my $numFiles=@files;
for(my $i=0 ; $i<$numFiles ; ++$i) {
  my $file=$files[$i];
  chomp $file;
  $file=~/([^\/]+)\.fastb/ || die $file;
  my $stem=$1;
  my $fastb=new Fastb("$inDir/$file");
  my $dnaTrack=$fastb->getTrackByName("DNA");
  my $transcript=$dnaTrack->getData();
  my $transcriptLen=length($transcript);
  my $sites=$sitesHash->{$stem};
  foreach my $site (@$sites) {
    my $begin=$site->getBegin();
    my $end=$site->getEnd();
    my $add=join(";",@{$site->{additionalFields}});
    $add=~/miRNA=([^;]+);/ || die $add;
    my $id=$1;
    $add=~/seq=([^;]+);/ || die $add;
    my $seed=$1;
    my $seedLen=length($seed);
    my $mature=$matureHash->{$id};
    my $matureLen=length($mature);
    my $rMature=Translation::reverseComplement(\$mature);

    my $seedPosInMature=findSeed($seed,$rMature);
    my $matureContextLen=$seedPosInMature;
    my $seedPosInTranscript=$begin;
    my $matureSeed=substr($rMature,$seedPosInMature,$seedLen);
    my $transcriptSeed=substr($transcript,$seedPosInTranscript,$seedLen);
    if($matureSeed ne $transcriptSeed && 
       substr($transcriptSeed,length($transcriptSeed)-1,1) ne "A")
      { die "$matureSeed vs. $transcriptSeed" }
    my $matureBeginInTranscript=$seedPosInTranscript-$matureContextLen;
    my $nickyBegin=$seedPosInTranscript-40;
    my $nickyEnd=$seedPosInTranscript+$seedLen+10;
    if($nickyBegin<0) { $nickyBegin=0 }
    if($nickyEnd>=$transcriptLen) { $nickyEnd=$transcriptLen-1 }
    my $nickyOffset=$seedPosInTranscript-$nickyBegin;
    my $nickySeq=substr($transcript,$nickyBegin,$nickyEnd-$nickyBegin+1);
    print ">$stem /miRNA=$id /seed=$seed /offset=$nickyOffset\n$nickySeq\n";
    #my $matureContext=substr();
    #my $transcriptContext=substr();

#===================================================================
goto SKIP;
    my $contextLen=22-$seedLen;
    my $twentyTwo=substr($transcript,$begin-$contextLen,22);
    my $shift=findShift($twentyTwo,$rMature);
    my $offset=22+$shift-$matureLen;
    if($offset>0) { $twentyTwo=substr($twentyTwo,$offset,22-$offset) }
    elsif($offset<0) { $rMature=substr($rMature,-$offset,$matureLen+$offset) }
    my $len22=length($twentyTwo);
    my $lenMature=length($rMature);
    my $contextLen=$len22-$seedLen;
    my $context1=substr($twentyTwo,0,$contextLen);
    my $context2=substr($rMature,0,$contextLen);
    my $sanity1=substr($twentyTwo,$contextLen,$seedLen);
    my $sanity2=substr($rMature,$contextLen,$seedLen);
    if(length($sanity1)!=length($sanity2)) {die "$twentyTwo $rMature"}
    if(substr($sanity1,0,length($sanity1)-1) ne substr($sanity2,0,length($sanity2)-1)) {die "$twentyTwo $rMature $sanity1 $sanity2"}
    my $matches=countMatches($context1,$context2);
    my $pctIdentity=$matches/$contextLen;
    my $alignment=$aligner->alignSeqs($context1,$context2,"DNA",1,1);
    my $alignIdent=$alignment->getPercentIdentity()/100;
    print "$pctIdentity\t$alignIdent\n";
#===================================================================
  SKIP:
  }
}

sub findSeed {
  my ($seed,$rna)=@_;
  if(substr($seed,length($seed)-1,1) eq "A") { chop $seed }
  my $L=length($rna);
  my $seedLen=length($seed);
  my $first=$L-8;
  my $last=$L-6;
  for(my $i=$first ; $i<=$last ; ++$i) {
    if($seed eq substr($rna,$i,$seedLen)) {return $i}
  }
  die;
}

sub countMatches {
  my ($A,$B)=@_;
  my $L=length($A);
  my $matches=0;
  for(my $i=0 ; $i<$L ; ++$i) {
    if(substr($A,$i,1) eq substr($B,$i,1)) { ++$matches }
  }
  return $matches;
}

sub alignTails {
  my ($A,$B,$aL,$bL,$shift)=@_;
  my $subLen=8-$shift;
  my $Bsub=substr($B,$bL-8,$subLen);
  my $Asub=substr($A,$aL-$subLen,$subLen);
  my $matches=0;
  for(my $i=0 ; $i<$subLen ; ++$i) {
    if(substr($Asub,$i,1) eq substr($Bsub,$i,1)) { ++$matches }
  }
  return $matches;
}

sub findShift {
  my ($A,$B)=@_;
  my $aL=length($A); my $bL=length($B);
  my $score1=alignTails($A,$B,$aL,$bL,0);
  my $score2=alignTails($A,$B,$aL,$bL,1);
  my $score3=alignTails($A,$B,$aL,$bL,2);
  my $shift;
  if($score1>$score2) {
    if($score1>$score3) {$shift=0}
    else {$shift=2}
  }
  elsif($score2>$score3) {$shift=1}
  else {$shift=2}
  return $shift;
}

