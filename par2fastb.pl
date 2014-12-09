#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use FastaReader;
use Translation;
use TempFilename;
use ParSite;
use ParGene;
use ParUTR;
use ParTranscript;
use ParGroup;

# THINGS YOU MIGHT WANT TO MODIFY:
my $MAX_GENES=-1; # set to -1 to make unlimited
my @keepSignalTypes=("Signal");#,"Background","ConversionPercent","ReadCount");
my $ANNOTATION_OUTFILE="paralyzer-peaks.gff";
my $INCLUDE_FOLD_SEQ=0;
my $MAX_CHROM_NAME_LEN=3;
my $MAX_GENE_LENGTH=1000000;
my $TWO_BIT_FILE="Genome.2bit";
my $DATASETS_SUBDIR="datasets";
my $NORMALIZATION_SUBDIR="normalization";
my $CLUSTERS_DIR="clusters";
my $NORM_PEAKS_TO_UNITY=1;
my $EXTEND_FOR_FOLDING=200; # number of bp to extend upstream
#my $UTRS_FILE="Hi.Expressed.3UTRs.txt";
my $UTRS_FILE="UTRs.txt";

#=====================================================================
# Don't modify anything below this point...
#=====================================================================
die "par2fastb.pl <outdir> <lib1> <lib2> ...\n" unless @ARGV>1;
my $BASEDIR=shift @ARGV;
my @libs=@ARGV;

my %keepSignalTypes;
foreach my $s (@keepSignalTypes) { $keepSignalTypes{$s}=1 }

# LOAD THE UTR FILE - - - - - - - - - - - - - - - - - - - - - - - - - -
print STDERR "loading utr coords\n";
my (%genes,%badGenes);
open(IN,$UTRS_FILE) || die;
while(<IN>) {
  chomp;
  my @fields=split/\t/,$_;
  my ($chr,$start,$end,$geneID,$strand,$transID)=@fields;
  next unless length($chr)<=$MAX_CHROM_NAME_LEN; ###
  my $gene=$genes{$geneID};
  if(!defined($gene))
    { $gene=$genes{$geneID}=new ParGene($chr,$strand) }
  if($chr ne $gene->getChr() || $strand ne $gene->getStrand())
    { $badGenes{$geneID}=1; next }
  my $utr=new ParUTR($start,$end);
  my $transcript=$gene->findOrCreateTranscript($transID);
  $transcript->addUTR($utr);
}
close(IN);
my @badGenes=keys %badGenes;
my $numBadGenes=@badGenes;
for(my $i=0 ; $i<$numBadGenes ; ++$i) { delete $genes{$badGenes[$i]} }

# SORT UTRS INTO TRANSLATION ORDER - - - - - - - - - - - - - - - - - - -
print STDERR "sorting utrs...\n";
my @geneIDs=keys %genes;
my $numGenes=@geneIDs;
for(my $i=0 ; $i<$numGenes ; ++$i) {
  my $geneID=$geneIDs[$i];
  my $gene=$genes{$geneID};
  my $strand=$gene->getStrand();
  my $transcriptIDs=$gene->getTranscriptIDs;
  foreach my $transID (@$transcriptIDs) {
    my $transcript=$gene->find($transID);
    $transcript->sortUTRs($strand);
  }
}

# LOAD THE CLUSTER FILES - - - - - - - - - - - - - - - - - - - - - - - -
print STDERR "loading clusters...\n";
foreach my $lib (@libs) {
  my $filename="$CLUSTERS_DIR/$lib.csv";
  if(!-e $filename) {print STDERR "WARNING: $filename not found\n"}
  next unless -e $filename;
  open(IN,$filename) || die "can't open file $filename\n";
  while(<IN>) {
    chomp;
    my @fields=split/,/,$_;
    my ($Chromosome,$Strand,$ClusterStart,$ClusterEnd,$ClusterID,$ClusterSequence,$ReadCount,$ModeLocation,$ModeScore,$ConversionLocationCount,$ConversionEventCount,$NonConversionEventCount,$FilterType,$TranscriptLocation,$TranscriptID,$GeneName)=@fields;
    next unless $TranscriptLocation eq "3UTR";
    $ModeScore=int(100*$ModeScore+5/9)/100;
    my $site=new ParSite($ClusterStart,$ClusterEnd,$ModeScore,$lib);
    my $gene=$genes{$GeneName}; ### NEVER FOUND
    next unless $gene;
    my $chr=$gene->getChr();
    my $strand=$gene->getStrand();
    if($Chromosome=~/chr(\S+)/) { $Chromosome=$1 }
    if($chr ne $Chromosome || $strand ne $Strand)
      {print STDERR "$chr vs $Chromosome, $strand vs $Strand, $GeneName\n";
       next }
    my $transcript=$gene->find($TranscriptID);
    next unless $transcript;
    my $numUTRs=$transcript->numUTRs();
    for(my $i=0 ; $i<$numUTRs ; ++$i) {
      my $utr=$transcript->getIthUTR($i);
      if($utr->containsInterval($ClusterStart,$ClusterEnd))
	{ $utr->addSite($site) }
    }
  }
  close(IN);
}

# LOAD THE LIBRARIES - - - - - - - - - - - - - - - - - - - - - - - - -
my %norm;
foreach my $libraryID (@libs) {
  print STDERR "processing library $libraryID\n";

  # LOAD GROUP-GENE MAPPING
  my %groupToGene;
  my $groupsFile="$DATASETS_SUBDIR/$libraryID.groups.csv";
  open(IN,$groupsFile) || die;
  while(<IN>) {
    chomp;
    my @fields=split/,/,$_;
    my $groupID=$fields[4];
    my $geneID=$fields[@fields-1];
    $groupToGene{$groupID}=$geneID;
  }
  close(IN);

  # Load normalization data
  my $infile="$NORMALIZATION_SUBDIR/$libraryID.norm";
  print STDERR "loading norm data $infile\n";
  open(IN,$infile) || die;
  while(<IN>) {
    if(/(chr\S+,\S,\d+,\d+),(\S+),(\S+)\s*$/) {
      my ($key,$signal,$val)=($1,$2,$3);
      $norm{$key}->{$signal}=$val;
    }
  }
  close(IN);

  # LOAD DISTRIBUTION FILE
  my $filename="$DATASETS_SUBDIR/$libraryID.distribution.csv";
  print STDERR "loading distr file $filename\n";
  open(IN,$filename) || die $filename;
  while(<IN>) {
    chomp;
    my @fields=split/,/,$_;
    my ($chr,$strand,$begin,$end,$groupID,$signalType)=@fields;
    next unless $keepSignalTypes{$signalType};
    my $geneID=$groupToGene{$groupID};
    my $gene=$genes{$geneID};
    next unless $gene;
    my $key="$chr,$strand,$begin,$end";
    die unless defined $norm{$key};###
    my $normConst=$norm{$key}->{$signalType};
    for(my $i=0 ; $i<6 ; ++$i) {shift @fields}
    my $data=[];
    @$data=@fields;
    my $n=@$data;
    if($normConst>0 && $normConst!~/Infinity/)
      { for(my $i=0 ; $i<$n ; ++$i) { $data->[$i]*=$normConst } }
    my $uniform=1;
    my $max=0;
    for(my $i=0 ; $i<$n ; ++$i) 
      { if($data->[$i]=~/nan/i) { $data->[$i]=0 } }
    for(my $i=0 ; $i<$n ; ++$i) { if($data->[$i]>$max) {$max=$data->[$i]} }
    for(my $i=0 ; $i<$n ; ++$i) { if($data->[$i]!=$max) {$uniform=0} }
    next if $uniform;
    if($max<=0) {die @$data}
    if($NORM_PEAKS_TO_UNITY && $max>0) 
      { for(my $i=0 ; $i<$n ; ++$i) { $data->[$i]/=$max } }
    my $group=new ParGroup($begin,$data,$libraryID,$signalType);
    $gene->addGroup($group);
  }
  close(IN);
}

removeEmptyTranscripts();

# GENERATE OUTPUT - - - - - - - - - - - - - - - - - - - - - - - - - - -
my $annotationOutfile=$ANNOTATION_OUTFILE;
if(@libs==1) {my $lib=$libs[0]; $annotationOutfile="$lib-peaks.gff"}
open(ANNO,">$annotationOutfile") || 
  die "can't write to file $ANNOTATION_OUTFILE";
my @geneIDs=keys %genes;
my $numGenes=@geneIDs;
for(my $i=0 ; $i<$numGenes ; ++$i) {
  if($MAX_GENES>=0 && $i>=$MAX_GENES) { last }
  if($i%100==0) {print STDERR "gene $i of $numGenes\n"}
  my $geneID=$geneIDs[$i];
  my $gene=$genes{$geneID};
  if($gene->isLongerThan($MAX_GENE_LENGTH)) {next}
  my $chr=$gene->getChr();
  my $strand=$gene->getStrand();
  my $groups=$gene->getGroups();
  if(!$groups) {die} ###
  my $transcriptIDs=$gene->getTranscriptIDs();
  my $numTranscripts=@$transcriptIDs;
  for(my $j=0 ; $j<$numTranscripts ; ++$j) {
    my $transcriptID=$transcriptIDs->[$j];
    my $transcript=$gene->find($transcriptID);
    my $numUTRs=$transcript->numUTRs();
    my %exonIDs;
    $transcript->sortUTRs($strand);
    for(my $k=0 ; $k<$numUTRs ; ++$k) {
      my $utr=$transcript->getIthUTR($k);
      if(!$utr) {die "FAIL TO DELETE"}
      my $utrBegin=$utr->getBegin()-1;###
      if(!$utrBegin) {die}
      my $utrEnd=$utr->getEnd();
      my $utrLen=$utrEnd-$utrBegin;
      my $isoID=$j;
      my $exonID=$k;
      my $outfile="$BASEDIR/$geneID\_$isoID\_$exonID.fastb";
      writeAnno(\*ANNO,$utr,$geneID,$isoID,$strand);
      my $tempFile=TempFilename::generate();
      system("rm $tempFile") if -e $tempFile;
      my $msg=`twoBitToFa $TWO_BIT_FILE $tempFile -seq=chr$chr -start=$utrBegin -end=$utrEnd -noMask`;
      if($msg=~/is not in/) {next}
      if(-z "$tempFile") {system("rm $tempFile"); next}
      my $reader=new FastaReader("$tempFile");
      my ($def,$seq)=$reader->nextSequence();
      if($strand eq "-") { $seq=Translation::reverseComplement(\$seq) }
      $reader->close();
      system("rm $tempFile");
      open(OUT,">$outfile") || die $outfile;
      #$seq="\U$seq";
      print OUT ">DNA /geneID=$geneID /transcriptID=$transcriptID /chr=$chr /begin=$utrBegin /end=$utrEnd /strand=$strand\n$seq\n";
      foreach my $libID (@libs) {
	foreach my $signalType (@keepSignalTypes) {
	  my @seq;
	  for(my $i=0 ; $i<$utrLen ; ++$i) {$seq[$i]=0}
	  foreach my $group (@$groups) {
	    next unless $group->getLibraryID() eq $libID &&
	      $group->getSignalType() eq $signalType;
	    my $groupBegin=$group->getBegin()-1;###
	    my $data=$group->getData();
	    my $dataLen=@$data;
	    my $pos=$groupBegin-$utrBegin;
	    for(my $i=0 ; $i<$dataLen && $pos<$utrLen ; ++$i, ++$pos) {
	      if($pos>=0) { $seq[$pos]=$data->[$i] }
	    }
	  }
	  print OUT "%$libID\-$signalType /length=$utrLen\n";
	  if($strand eq "-") {
	    my @tmp;
	    for(my $i=0 ; $i<$utrLen ; ++$i) {unshift @tmp,$seq[$i]}
	    @seq=@tmp;
	  }
	  for(my $i=0 ; $i<$utrLen ; ++$i) {print OUT "$seq[$i]\n"}
	}
      }
      close(OUT);
    }
  }
}
close(ANNO);

#########################################################################
sub writeAnno
  {
    my ($fh,$utr,$geneID,$isoID,$strand)=@_;
    my $utrBegin=$utr->getBegin();
    my $utrEnd=$utr->getEnd();
    my $utrOffset=$utr->getTranscriptOffset();
    my $numSites=$utr->getNumSites();
    for(my $i=0 ; $i<$numSites ; ++$i) {
      my $site=$utr->getIthSite($i);
      my $score=$site->getScore();
      my $lib=$site->getLibraryID();
      my $length=$site->getLength();
      my $begin=$strand eq "+" ?
	$site->getBegin()-$utrBegin+$utrOffset :
	  $utrEnd-$site->getEnd()+$utrOffset;
      ++$begin;
      my $end=$begin+$length-1;

      # GFF is one-based & inclusive:
      print $fh "$geneID\_$isoID\tparalyzer\t$lib\t$begin\t$end\t$score\t$strand\t.\n";
    }
  }

sub removeEmptyTranscripts
  {
    my @geneIDs=keys %genes;
    my $numGenes=@geneIDs;
    for(my $i=0 ; $i<$numGenes ; ++$i) {
      my $geneID=$geneIDs[$i];
      my $gene=$genes{$geneID};
      my $groups=$gene->getGroups();
      if(!$groups) { delete $genes{$geneID}; next }
      my $transcriptIDs=$gene->getTranscriptIDs();
      foreach my $transcriptID (@$transcriptIDs) {
	my $transcript=$gene->find($transcriptID);
	my $transcriptOK=0;
	foreach my $group (@$groups) {
	  if($transcript->overlapsGroup($group))
	    { $transcriptOK=1 ; last }
	}
	if(!$transcriptOK) { $gene->removeTranscript($transcriptID) }
      }
    }
  }

