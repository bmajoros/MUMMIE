#!/usr/bin/perl
use strict;
use FastaReader;
use ProgramName;
$|=1;

my $name=ProgramName::get();
die "$name <var> <VITERBI/POSTERIOR> <shuffled-seeds.fasta> <unshuffled-seeds.fasta> <expression-levels.txt> <FASTB-dir> <foreground-states> <metamodel.txt> <submodels.txt> <schema.txt> <groups.gff> <bg.hmm> <signal-name> <output-dir>\n"
  unless @ARGV==12;
my ($VARIANCE,$decoder,$SHUFFLED,$UNSHUFFLED,$EXPRESSION_LEVELS,
    $TEST_DIR,$STATES,$METAMODEL,$SUBMODELS,$SCHEMA,$GROUPS,$BGHMM,$SIGNAL,
    $OUTDIR)=@ARGV;

my $MUMMIE="/home/ohler/bmajoros/MUMMIE";
my $SKIPTRANS=0.0;
my $USE_CLUSTERS=0;
my $DASH_G=($decoder eq "VITERBI" ? "-g" : "-I");
my $EVAL_SCRIPT=($decoder eq "VITERBI" ? "eval-SnSp.pl" : "eval-SnSp2.pl");
my $NUM_ITERATIONS=100;
my $TRAINING_SEEDS="training-seeds.fasta";
my $PEAK_THRESHOLD=0.999;
my $TRANS_6MER=0.33;#0.05;
my $TRANS_7MER=0.33;#0.19;
my $TRANS_8MER=0.33;#0.76;

#if(-e "$OUTDIR") { system("rm -r $OUTDIR") }
system("mkdir -p $OUTDIR");

print "finding peaks...\n";
system("$MUMMIE/find-peaks -t $PEAK_THRESHOLD -d $TEST_DIR $SCHEMA $SIGNAL > $OUTDIR/WT-peaks.gff");

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
  $def=~/>(\S+)_sh/ || die;
  my $id=$1;
  $seq=~s/U/T/g;
  push @{$shuffled{$id}},$seq;
}
$reader->close();

my $TP=0;
my $FP=0;
my $naiveTP=0;
my $naiveFP=0;

print "writing .q files...\n";
for(my $i=0 ; $i<$NUM_ITERATIONS ; ++$i) {
  my $iter=$i;
  my $X="$OUTDIR/$i";

  # First, sample a set of 35 shuffled seeds, one per unshuffled seed
  my @sample;
  foreach my $pair (@unshuffled) {
    push @sample,$pair;
    my ($id,$seq)=@$pair;
    my $array=$shuffled{$id};
    my $n=@$array;
    my $index=int(rand($n));
    my $shuff=$array->[$index];
    push @sample,["$id\_sh_$index",$shuff];
  }

  # Write out the training set
  open(OUT,">$X.$TRAINING_SEEDS") || die;
  foreach my $pair (@sample) {
    my ($def,$seq)=@$pair;
    print OUT ">$def\n$seq\n";
  }
  close(OUT);

  open(OUT,">$X.q") || die;
  my $clusters=($USE_CLUSTERS ? "WT.clusters.gff" : "");
  print OUT "#!/bin/tcsh
#
#\$ -S /bin/tcsh -cwd
#\$ -o $OUTDIR/$iter.mpirun -j y
#\$ -l highprio
#\$ -N shf$iter
#\$ -q *\@ohler*,*\@igspnih-n1*,*\@igspnih-n3*,*\@igspnih-n4*,*\@igspnih-n6*,*\@igspncbc*
cd $OUTDIR
hostname
$MUMMIE/motif-scan.pl $X.$TRAINING_SEEDS $TEST_DIR 0 > $X.scan.gff
echo trace1
$MUMMIE/nicky/sites-closest-to-a-peak.pl $OUTDIR/WT-peaks.gff $X.scan.gff $clusters > $X.proximal.gff
echo trace2
$MUMMIE/nicky/eval-SnSp.pl $X.proximal.gff $GROUPS > $X.naive.eval
echo trace3
$MUMMIE/nicky/make-42state-site.pl $X.$TRAINING_SEEDS $TRANS_6MER $TRANS_7MER $TRANS_8MER 4 $SCHEMA 7 $BGHMM $X.site.hmm $EXPRESSION_LEVELS
echo trace4
cp $SUBMODELS $X.submodels
echo trace5
/home/ohler/bmajoros/bin/sub $X.submodels site.hmm $X.site.hmm
echo trace6
$MUMMIE/model-combiner $METAMODEL $X.submodels $X.hmm
$MUMMIE/hmm-edit $X.hmm DTRK phastcons
$MUMMIE/hmm-edit $X.hmm VAR all 0 $VARIANCE
echo trace7
rm $X.submodels
echo trace8
$MUMMIE/parse $X.hmm $TEST_DIR -d -p $DASH_G $STATES > $X.predictions.gff.tmp
echo trace9
$MUMMIE/nicky/identify-miRNAs.pl $X.predictions.gff.tmp $X.$TRAINING_SEEDS > $X.predictions.gff.tmp2
rm $X.predictions.gff.tmp
$MUMMIE/combine-miRNA-predictions.pl $X.predictions.gff.tmp2 > $X.predictions.gff
rm $X.predictions.gff.tmp2
echo trace10
$EVAL_SCRIPT $X.predictions.gff $GROUPS > $X.mummie.eval
echo done
";
  close(OUT);
  system("qsub $X.q >& /dev/null");
}

print "\nall jobs have been submitted\n";

