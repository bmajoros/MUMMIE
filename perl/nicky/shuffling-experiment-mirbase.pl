#!/usr/bin/perl
use strict;
use FastaReader;
use ProgramName;
$|=1;

my $SKIPTRANS=0.0;
my $USE_CLUSTERS=0;
my $VARIANCE=0.25; # 0.04 is the "standard" setting for model N3
my $DASH_G="-g"; # "-g" or "-I";
my $EVAL_SCRIPT="eval-SnSp.pl"; # "eval-SnSp2.pl";
my $EXPRESSION_LEVELS="\"\"";
#my $EXPRESSION_LEVELS="expression-levels.txt";

my $NUM_ITERATIONS=100;
#my $SHUFFLED="nicky-top50-shuffled-seeds-bygroup.fasta";
#my $UNSHUFFLED="nicky-top50-seeds.fasta";
my $TRAINING_SEEDS="training-seeds.fasta";
my $METAMODEL="shuffled-metamodel.txt"; # N3
my $SUBMODELS="shuffled-submodels.txt"; # N3
my $TEST_DIR="fastb-WT-nocons-longest";
my $SCHEMA="WT-cons.schema";
my $STATES="5-45"; # N3 w/42-state site
my $PEAK_THRESHOLD=0.999;
my $TRANS_6MER=0.33;#0.05;
my $TRANS_7MER=0.33;#0.19;
my $TRANS_8MER=0.33;#0.76;

my $name=ProgramName::get();
die "$name <seeds.fasta> <shuffled-seeds.fasta> <basedir>\n" unless @ARGV==3;
my ($UNSHUFFLED,$SHUFFLED,$BASEDIR)=@ARGV;

if(-e "$BASEDIR") { system("rm -r $BASEDIR") }
system("mkdir $BASEDIR");

print "finding peaks...\n";
system("MUMMIE/find-peaks -t $PEAK_THRESHOLD -d $TEST_DIR $SCHEMA AGO-Signal > WT-peaks.gff");

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
  my $X="$BASEDIR/$i";

  # First, sample a set of 35 shuffled seeds, one per unshuffled seed
  my @sample;
  foreach my $pair (@unshuffled) {
    push @sample,$pair;
    my ($id,$seq)=@$pair;
    my $array=$shuffled{$id};
    die "undefined: $id" unless $array;
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
#\$ -o /home/ohler/bmajoros/MUMMIE/nicky/$X.mpirun -j y
#\$ -l highprio
#\$ -N shuf$iter
#\$ -q *\@ohler*,*\@igspnih-n1*,*\@igspnih-n3*,*\@igspnih-n4*,*\@igspnih-n6*,*\@igspncbc*
cd /home/ohler/bmajoros/MUMMIE/nicky
hostname
MUMMIE/motif-scan.pl $X.$TRAINING_SEEDS $TEST_DIR 0 > $X.scan.gff
echo trace1
sites-closest-to-a-peak.pl WT-peaks.gff $X.scan.gff $clusters > $X.proximal.gff
echo trace2
eval-SnSp.pl $X.proximal.gff groups.gff > $X.naive.eval
echo trace3
make-42state-site.pl $X.$TRAINING_SEEDS $TRANS_6MER $TRANS_7MER $TRANS_8MER 4 WT-cons.schema 7 bg.hmm $X.site.hmm $EXPRESSION_LEVELS
echo trace4
cp $SUBMODELS $iter.submodels
echo trace5
sub $iter.submodels site.hmm $X.site.hmm
echo trace6
MUMMIE/model-combiner $METAMODEL $iter.submodels $X.hmm
MUMMIE/hmm-edit $X.hmm DTRK phastcons
MUMMIE/hmm-edit $X.hmm VAR all 0 $VARIANCE
cho trace7
rm $iter.submodels
echo trace8
MUMMIE/parse $X.hmm $TEST_DIR -d -p $DASH_G $STATES > $X.predictions.gff.tmp
echo trace9
identify-miRNAs.pl $X.predictions.gff.tmp $X.$TRAINING_SEEDS > $X.predictions.gff.tmp2
rm $X.predictions.gff.tmp

MUMMIE/combine-miRNA-predictions.pl $X.predictions.gff.tmp2 > $X.predictions.gff
#cp $X.predictions.gff.tmp2 $X.predictions.gff

rm $X.predictions.gff.tmp2
echo trace10
# eval-experiment2.pl $X.predictions.gff > $X.mummie.eval
$EVAL_SCRIPT $X.predictions.gff groups.gff > $X.mummie.eval
echo trace11
";
  close(OUT);

  system("qsub $X.q >& /dev/null");
}


print "\nall jobs have been submitted\n";

