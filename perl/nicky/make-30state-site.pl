#!/usr/bin/perl
use strict;
use ProgramName;
use FastaReader;
use TempFilename;

my $name=ProgramName::get();
die "$name <8mers.fasta> <6mer-prob> <7mer-prob> <8mer-prob> <num-mix-comp> <schema> <order> <out.hmm>\n" unless @ARGV==8;
my ($seedsFile,$sixProb,$sevenProb,$eightProb,$mix,$schema,$order,$outfile)=
  @ARGV;

my $SUBDIR=TempFilename::generate();
system("mkdir $SUBDIR");

# Load 8mers
my @seeds;
my $reader=new FastaReader($seedsFile);
while(1) {
  my ($def,$seq)=$reader->nextSequence();
  last unless $def;
  push @seeds,[$def,$seq];
}

# Make canonical models for each length
my $nextSeed=int(rand(1000000));
for(my $L=6 ; $L<=8 ; ++$L) {
  my $numStates=$L+1;
  my $outfile="$SUBDIR/canonical$L.hmm";
  system("MUMMIE/random-HMM -c V -u -s $nextSeed $numStates 1 $mix $schema $order $outfile");
  die "MUMMIE/random-HMM -c V -u -s $nextSeed $numStates 1 $mix $schema $order $outfile" unless -e $outfile;
  ++$nextSeed;
}
my $tgfFile="$SUBDIR/tgf.tgf";
system("MUMMIE/make-tgf.pl $schema full $tgfFile");

# Train each submodel
for(my $L=6 ; $L<=8 ; ++$L) {
  my $inHMM="$SUBDIR/canonical$L.hmm";
  my $outHMM="$SUBDIR/canonical$L-trained.hmm";
  my @offsets=(0);
  if($L<8) { push @offsets,1 }
  my $nmersFile=makeTrainingFile(\@seeds,\@offsets,$L);
  my $nmersDir="$SUBDIR/nmers$L";#TempFilename::generate();
  my $tieProfile="$SUBDIR/tie$L.profile";#TempFilename::generate();
  makeTieProfile($tieProfile,$L);
  system("mkdir $nmersDir");
  system("MUMMIE/fasta-to-fastb.pl $nmersFile $schema $nmersDir");
  system("MUMMIE/baum-welch $inHMM $tgfFile $nmersDir 2 $outHMM -R -s $nextSeed -u -t $tieProfile -n 0.00001");
  system("rm $nmersFile $tieProfile");
  system("rm -r $nmersDir");
  ++$nextSeed;
}
my $inHMM="$SUBDIR/canonical7.hmm";
my $outHMM="$SUBDIR/canonical7-trained2.hmm";
my $nmersFile=makeTrainingFile(\@seeds,[1],7);
my $nmersDir="$SUBDIR/nmers";#TempFilename::generate();
my $tieProfile="$SUBDIR/tie.profile";#TempFilename::generate();
makeTieProfile($tieProfile,7);
system("mkdir $nmersDir");
system("MUMMIE/fasta-to-fastb.pl $nmersFile $schema $nmersDir");
system("MUMMIE/baum-welch $inHMM $tgfFile $nmersDir 2 $outHMM -R -s $nextSeed -u -t $tieProfile -n 0.00001");
system("rm $nmersFile $tieProfile");
system("rm -r $nmersDir");
++$nextSeed;

# Combine submodels into one big model
system("cd $SUBDIR ; ~/MUMMIE/hmm-extract-state canonical6-trained.hmm all 6mer");
system("cd $SUBDIR ; ~/MUMMIE/hmm-extract-state canonical7-trained.hmm all 7mer");
system("cd $SUBDIR ; ~/MUMMIE/hmm-extract-state canonical7-trained2.hmm all 7mer2");
system("cd $SUBDIR ; ~/MUMMIE/hmm-extract-state canonical8-trained.hmm all 8mer");
open(OUT,">$SUBDIR/metamodel30.txt") || die;
my $halfSeven=$sevenProb/2;
print OUT
"0 -> 1 : $sixProb
0 -> 7 : $halfSeven
0 -> 14 : $halfSeven
0 -> 22 : $eightProb
";
linkStates(\*OUT,1,6);
linkStates(\*OUT,7,13);
linkStates(\*OUT,14,19);
linkStates(\*OUT,22,28);
print OUT
"6 -> 0 : 1
13 -> 0 : 1
19 -> 20 : 0.5
19 -> 21 : 0.5
20 -> 0 : 1
21 -> 0 : 1
28 -> 21 : 0.5
28 -> 29 : 0.5
29 -> 0 : 1
";
close(OUT);
makeAfiles("$SUBDIR/A30.hmm",$schema);
open(OUT,">$SUBDIR/submodels30.txt") || die;
addSubmodels(\*OUT,"6mer",1,6);
addSubmodels(\*OUT,"7mer",7,13);
addSubmodels(\*OUT,"7mer2",14,20);
print OUT "21 = A30.hmm\n";
addSubmodels(\*OUT,"8mer",22,29);
close(OUT);
system("rm $tgfFile");
my $pwd=`pwd`;
chomp $pwd;
system("cd $SUBDIR ; ~/MUMMIE/model-combiner metamodel30.txt submodels30.txt $pwd/$outfile");

# Set mixture parameters
system("MUMMIE/hmm-edit $outfile MIX ALL 0 0 MIX ALL 1 0 MIX ALL 2 0 MIX ALL 3 1");
system("MUMMIE/hmm-edit $outfile MEAN 0 0 1e-7   MEAN 0 1 0.531313");
system("MUMMIE/hmm-edit $outfile MEAN 1 0 0.5    MEAN 1 1 0.531313");
system("MUMMIE/hmm-edit $outfile MEAN 2 0 1.0    MEAN 2 1 0.531313");
system("MUMMIE/hmm-edit $outfile MEAN 3 0 0.5    MEAN 3 1 0.993844");
system("MUMMIE/hmm-edit $outfile COV ALL ALL ALL 0");
system("MUMMIE/hmm-edit $outfile VAR  0 0 0.01   VAR  0 1 1.10313");
system("MUMMIE/hmm-edit $outfile VAR  1 0 0.01   VAR  1 1 1.10313");
system("MUMMIE/hmm-edit $outfile VAR  2 0 0.01   VAR  2 1 1.10313");
system("MUMMIE/hmm-edit $outfile VAR  3 0 0.01   VAR  3 1 1.15067");

# Clean up
system("rm -r $SUBDIR");

#=====================================================================
sub makeAfiles {
  my ($outHMM,$schema)=@_;
  my $fasta=TempFilename::generate();
  my $dir=TempFilename::generate();
  system("mkdir $dir");
  open(OUT,">$fasta") || die;
  print OUT ">1\nA\n>2\nA\n>3\nA\n";
  close(OUT);
  system("MUMMIE/fasta-to-fastb.pl $fasta $schema $dir");
  my $inHMM=TempFilename::generate();
  system("MUMMIE/random-HMM -c 0 2 1 $mix $schema $order $inHMM");
  system("MUMMIE/hmm-edit $inHMM TRANS 1 1 0");
  my $tieProfile=TempFilename::generate();
  makeTieProfile($tieProfile,1);
  system("MUMMIE/baum-welch $inHMM $tgfFile $dir 2 $outHMM -R -s $nextSeed -u -t $tieProfile -n 0.00001");
  ++$nextSeed;
  system("rm $fasta $inHMM $tieProfile ; rm -r $dir");
}
sub addSubmodels {
  my ($fh,$tag,$begin,$end)=@_;
  my $L=$end-$begin+1;
  for(my $q=1 ; $q<=$L ; ++$q) {
    my $r=$q+$begin-1;
    print $fh "$r = $tag$q.hmm\n";
  }
}
sub linkStates {
  my ($fh,$begin,$end)=@_;
  for(my $q=$begin ; $q<$end ; ++$q) {
    my $nextQ=$q+1;
    print $fh "$q -> $nextQ : 1\n";
  }
}
sub makeTrainingFile {
  my ($seeds,$offsets,$L)=@_;
  my $filename=TempFilename::generate();
  open(OUT,">$filename") || die;
  foreach my $pair (@$seeds) {
    my ($def,$seed)=@$pair;
    $def=~/>(\S+)/;
    my $id=$1;
    foreach my $offset (@$offsets) {
      my $subseed=substr($seed,$offset,$L);
      my $def=">$id:$offset-$L";
      print OUT "$def\n$subseed\n";
    }
  }
  close(OUT);
  return $filename;
}
sub makeTieProfile {
  my ($filename,$L)=@_;
  open(OUT,">$filename") || die $filename;
  print OUT "fix means\nfix covariance_matrix\nfix transitions\n";
  print OUT "fix weights in states ";
  for(my $q=1 ; $q<$L ; ++$q) { print OUT "$q," }
  print OUT "$L\n";
  close(OUT);
}




