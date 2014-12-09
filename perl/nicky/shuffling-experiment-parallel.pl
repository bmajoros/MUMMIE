#!/usr/bin/perl
use strict;
use FastaReader;
use ProgramName;
use SGE;
$|=1;
my $name=ProgramName::get();
die "$name <var> <skip-trans-P>\n" unless @ARGV==2;
my ($VARIANCE,$SKIPTRANS)=@ARGV;
#my $POSTERIOR_THRESHOLD=$SKIPTRANS; # hack

# SCORE OPTIONS ARE: score3p, scoreAU, scorePos, scoreAbund,
#                    scoreSPS, branchLen, score
#my $PARALLEL="fc/shuffling/parallel-PD$POSTERIOR_THRESHOLD";
my $PARALLEL="fc/shuffling/parallel-var$VARIANCE";
my $NUM_ITERATIONS=1000;
my $RUN_MUMMIE=1; # 0=no, 1=yes
my $MAX_JOBS=2000;
my $TARGETSCAN_SCORETYPE="branchLen";
my $TARGETSCAN_VARIANCE=3.0;
my $TARGETSCAN_BG_MEAN=0.0;
my $TARGETSCAN_FG_MEAN=3.0;
my $NUM_SEEDS=50;
my $TRANS_6MER=0.33;#0.05;
my $TRANS_7MER=0.33;#0.75;#0.19;
my $TRANS_8MER=0.33;#0.20;#0.76;
#my $EXPRESSION_LEVELS="expression-levels.txt";
my $EXPRESSION_LEVELS="\"\"";
my $USE_CLUSTERS=1; # only for the naive scan
my $TS_SEED_TYPES="N"; # Y or N
my $UNIFORM_SEEDS="N";
my $SOURCE_OF_BLS_SCORES="15percent/targetscan/all.gff";
#my $SOURCE_OF_BLS_SCORES="fc/datasets/modified-BLS/all.gff";

system("MUMMIE/get-groups.pl AGO-Signal fastb-WT-nocons-longest 0 > groups.gff");

#==============
my $UNSHUFFLED="15percent/mirbase/top$NUM_SEEDS.fasta";
my $SHUFFLED="15percent/shuffles/top$NUM_SEEDS.fasta";
my $TEST_DIR="fastb-WT-nocons-longest";
#my $TEST_DIR="fastb-WT-targetscan";
my $SCHEMA="WT-cons.schema";
#my $SCHEMA="targetscan.schema";

my $DASH_G="-g"; # "-g" or "-I";
my $EVAL_SCRIPT="eval-SnSp.pl"; # "eval-SnSp2.pl";
my $TRAINING_SEEDS="training-seeds.fasta";
my $METAMODEL="shuffled-metamodel.txt"; # N3
my $SUBMODELS="shuffled-submodels.txt"; # N3
my $STATES="5-45"; # N3 w/42-state site
my $PEAK_THRESHOLD=0.999;
my $sge=new SGE();

if(-e "$PARALLEL") { system("rm -r $PARALLEL") }
system("mkdir -p $PARALLEL");

system("MUMMIE/find-peaks -t $PEAK_THRESHOLD -d $TEST_DIR $SCHEMA AGO-Signal > $PARALLEL/WT-peaks.gff");

# Load the unshuffled seeds
my @unshuffled;
my $reader=new FastaReader($UNSHUFFLED);
while(1) {
  my ($def,$seq)=$reader->nextSequence();
  last unless $def;
  $def=~/>(\S+)/ || die;
  my $id=$1;
  push @unshuffled,[$id,$seq];
}
$reader->close();

# Load all shuffled seeds
my %shuffled;
my $reader=new FastaReader($SHUFFLED);
while(1) {
  my ($def,$seq)=$reader->nextSequence();
  last unless $def;
  $def=~/>(\S+)/ || die;
  my $originalID=$1;
  $def=~/>(\S+)_sh/ || die;
  my $id=$1;
  $seq=~s/U/T/g;
  push @{$shuffled{$id}},[$originalID,$seq];
}
$reader->close();

my $TP=0;
my $FP=0;
my $naiveTP=0;
my $naiveFP=0;

for(my $i=0 ; $i<$NUM_ITERATIONS ; ++$i) {
  my $iter=$i;
  my $X="$PARALLEL/$i";

  # First, sample a set of shuffled seeds, one per unshuffled seed
  my @sample;
  foreach my $pair (@unshuffled) {
    push @sample,$pair;
    my ($id,$seq)=@$pair;
    my $array=$shuffled{$id};
    if(!$array) {die "no shuffles for $id"}
    my $n=@$array;
    my $index=int(rand($n));
    my $shuff=$array->[$index];
    push @sample,$shuff;
    #push @sample,["$id\_sh_$index",$shuff];
  }

  # Write out the training set
  open(OUT,">$X.$TRAINING_SEEDS") || die;
  foreach my $pair (@sample) {
    my ($def,$seq)=@$pair;
    print OUT ">$def\n$seq\n";
  }
  close(OUT);

  next unless $RUN_MUMMIE;
  open(OUT,">$X.q") || die;
  my $clusters=($USE_CLUSTERS ? "WT.clusters.gff" : "");
  print OUT "#!/bin/tcsh
#
#\$ -S /bin/tcsh -cwd
#\$ -o /home/ohler/bmajoros/MUMMIE/nicky/$X.mpirun -j y
#\$ -l highprio
#\$ -N shufex$iter
#\$ -l scr_free=1G
#\$ -q *\@ohler*,*\@igspnih-n1*,*\@igspnih-n3*,*\@igspnih-n4*,*\@igspnih-n6*,*\@igspncbc*
cd /home/ohler/bmajoros/MUMMIE/nicky
hostname
MUMMIE/motif-scan.pl $X.$TRAINING_SEEDS $TEST_DIR 0 > $X.scan.gff
echo trace1
sites-closest-to-a-peak.pl $PARALLEL/WT-peaks.gff $X.scan.gff $clusters > $X.proximal.gff
rm $X.scan.gff
echo trace2
eval-SnSp.pl $X.proximal.gff groups.gff > $X.naive.eval
echo trace3
rm $X.proximal.gff

exit

make-42state-site.pl $X.$TRAINING_SEEDS $TRANS_6MER $TRANS_7MER $TRANS_8MER 4 WT-cons.schema 7 bg.hmm $X.site.hmm $EXPRESSION_LEVELS $TS_SEED_TYPES $UNIFORM_SEEDS
echo trace4
cp $SUBMODELS $X.submodels
echo trace5
sub $X.submodels site.hmm $X.site.hmm
echo trace6
MUMMIE/model-combiner $METAMODEL $X.submodels $X.hmm
rm $X.submodels
rm $X.site.hmm
MUMMIE/hmm-edit $X.hmm TRANS 4 46 $SKIPTRANS
MUMMIE/hmm-edit $X.hmm DTRK phastcons
MUMMIE/hmm-edit $X.hmm VAR all 0 $VARIANCE
echo trace7

### TARGETSCAN PART (CONSERVATION: BRANCH LENGTH SCORE)
#MUMMIE/hmm-edit $X.hmm TRK targetscan
#echo trace7a
#MUMMIE/hmm-edit $X.hmm MEAN all -- 1 $TARGETSCAN_BG_MEAN
#echo trace7b
#MUMMIE/hmm-edit $X.hmm VAR  all -- 1 $TARGETSCAN_VARIANCE
#echo trace7c
#MUMMIE/hmm-edit $X.hmm MEAN 3   -- 1 $TARGETSCAN_FG_MEAN
#echo trace7d
#if( -e /scratch/bmajoros/$iter.chunks ) then
#  rm -r /scratch/bmajoros/$iter.chunks
#endif
#mkdir -p /scratch/bmajoros/$iter.chunks
#filter-gff-by-miRNAs.pl $SOURCE_OF_BLS_SCORES $X.training-seeds.fasta | sort > $X.targetscan-filtered.gff
#echo trace7f
#add-targetscan-tracks.pl $X.targetscan-filtered.gff fastb-WT-nocons-longest /scratch/bmajoros/$iter.chunks $TARGETSCAN_SCORETYPE
#rm $X.targetscan-filtered.gff
####

echo trace8
MUMMIE/parse $X.hmm $TEST_DIR -d -p $DASH_G $STATES > $X.predictions.gff.tmp
#MUMMIE/parse $X.hmm /scratch/bmajoros/$iter.chunks -d -p $DASH_G $STATES > $X.predictions.gff.tmp

echo trace9
identify-miRNAs.pl $X.predictions.gff.tmp $X.$TRAINING_SEEDS > $X.predictions.gff.tmp2
rm $X.predictions.gff.tmp
MUMMIE/combine-miRNA-predictions.pl $X.predictions.gff.tmp2 > $X.predictions.gff
#gff-threshold-scores.pl $X.predictions.gff \$POSTERIOR_THRESHOLD > $X.predictions.gff.tmp2
#mv $X.predictions.gff.tmp2 $X.predictions.gff
rm $X.predictions.gff.tmp2
echo trace10
$EVAL_SCRIPT $X.predictions.gff groups.gff > $X.mummie.eval
# rm $X.hmm $X.*hmm $X.q $X.mpirun $X.*gff
echo done
";
  close(OUT);

  system("qsub $X.q >& /dev/null");
  while($sge->countJobs()>=$MAX_JOBS) { sleep(15) }
}


print "all jobs have been submitted for variance $VARIANCE\n";

