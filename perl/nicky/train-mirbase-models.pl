#!/usr/bin/perl
use strict;
use ProgramName;

my @files=`ls mirbase/unshuffled`;
foreach my $file (@files) {
  next unless $file=~/(\S+)\.txt/;
  my $stem=$1;
  my $shuffStem="${stem}sh";
  my $rawUnshuffled="mirbase/unshuffled/$stem.txt";
  my $rawShuffled="mirbase/group/$shuffStem.txt";
  my $unshuffledSeeds="mirbase/$stem-seeds.fasta";
  my $shuffledSeeds="mirbase/$stem-shuffled-seeds.fasta";
  my $unshuffledHMM="mirbase/$stem.hmm";
  my $unshuffledPredictions="mirbase/$stem-Viterbi.gff";

goto DECODE;
  system("cat $rawUnshuffled | awk '{print \">\" \$1 \"\\n\" \$2}' > mirbase/mature.fasta");
  system("MUMMIE/get-mirna-training-set.pl mirbase/mature.fasta 0 dummy.file $unshuffledSeeds");
  system("cat $rawShuffled | awk '{print \">\" \$1 \"\\n\" \$2}' > mirbase/mature.fasta");
  system("MUMMIE/get-mirna-training-set.pl mirbase/mature.fasta 0 dummy.file $shuffledSeeds");
 DECODE:

goto SKIP;
  System("make-42state-site.pl $unshuffledSeeds 0.33 0.33 0.33 4 WT-cons.schema 7 bg.hmm site42.hmm \"\"");
  System("MUMMIE/model-combiner metamodel42.txt submodels42.txt $unshuffledHMM");
  System("MUMMIE/hmm-edit $unshuffledHMM DTRK phastcons");

#goto SKIP;
  System("MUMMIE/parse $unshuffledHMM fastb-WT-nocons-longest -p -g 5-45 -d > $unshuffledPredictions.tmp");
  System("identify-miRNAs.pl $unshuffledPredictions.tmp $unshuffledSeeds > $unshuffledPredictions ; rm $unshuffledPredictions.tmp");
  System("MUMMIE/parse $unshuffledHMM fastb-WT-nocons-longest -p -I 5-45 -d > $unshuffledPredictions.tmp");
  System("identify-miRNAs.pl $unshuffledPredictions.tmp $unshuffledSeeds > mirbase/$stem-PD.gff ; rm $unshuffledPredictions.tmp");
  next;

 SKIP:
  my $baseDir="mirbase/shuffling-experiments/$stem";
  system("shuffling-experiment-mirbase.pl $unshuffledSeeds $shuffledSeeds $baseDir");
}

sub System {
  my ($cmd)=@_;
  print "$cmd\n";
  system($cmd);
}

