#!/usr/bin/perl
use strict;
use FastaReader;
$|=1;

my $SKIPTRANS=0.0;
my $USE_CLUSTERS=1;
my $VARIANCE=0.1; # 0.04 is the "standard" setting for model N3
my $DASH_G="-I"; # "-g";
my $EVAL_SCRIPT="eval-SnSp2.pl"; # "eval-SnSp.pl";
my $NUM_ITERATIONS=10;

#my $SHUFFLED="shuffled35-unique.fasta";
#my $UNSHUFFLED="top35-unique.fasta";
#my $TRAINING_SEEDS="training-seeds.fasta";
my $SHUFFLED="nicky-top50-shuffled-seeds.fasta";
my $UNSHUFFLED="nicky-top50-seeds.fasta";
my $TRAINING_SEEDS="training-seeds.fasta";
my $METAMODEL="shuffled-metamodel.txt"; # N3
my $SUBMODELS="shuffled-submodels.txt"; # N3
#my $METAMODEL="N5-metamodel.txt";
#my $SUBMODELS="N5-shuffled-submodels.txt";
#my $SUBMODELS="N5c-submodels-shuffled.txt";
#my $METAMODEL="N4-metamodel.txt";
#my $SUBMODELS="N4-shuffled-submodels.txt";
#my $TEST_DIR="fastb-WT-cons-longest";
#my $TEST_DIR="fastb-WT-fold-longest";
#my $TEST_DIR="fastb-WT-primates-longest";
my $TEST_DIR="fastb-WT-nocons-longest";
#my $TEST_DIR="fastb-WT-fold-longest";
my $SCHEMA="WT-cons.schema";
#my $SCHEMA="WT-fold.schema";
#my $SCHEMA="WT-fold.schema";
#my $STATES="4-32,37-65"; # N5
#my $STATES="5-33"; # N3
#my $STATES="5-13,15-28,31-36,38-43"; # N3 w/42-state site
my $STATES="5-45"; # N3 w/42-state site
#my $STATES="3-31"; # N4
my $PEAK_THRESHOLD=0.999;
my $TRANS_6MER=0.33;#0.05;
my $TRANS_7MER=0.33;#0.19;
my $TRANS_8MER=0.33;#0.76;

if(-e "parallel") { system("rm -r parallel") }
system("mkdir parallel");

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
  my $X="parallel/$i";

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

  my $sampleIndex=0;
  foreach my $pair (@sample) {
    my ($def,$seq)=@$pair;
    ++$sampleIndex;
    my $X="parallel/$iter.$sampleIndex";
    open(OUT,">$X.$TRAINING_SEEDS") || die;
    print OUT ">$def\n$seq\n";
    close(OUT);
    open(OUT,">$X.q") || die;
    my $clusters=($USE_CLUSTERS ? "WT.clusters.gff" : "");
    print OUT "#!/bin/tcsh
#
#\$ -S /bin/tcsh -cwd
#\$ -o /home/ohler/bmajoros/MUMMIE/nicky/$X.mpirun -j y
#\$ -l highprio
#\$ -N s$iter\_$sampleIndex
#\$ -q *\@ohler*,*\@igspnih-n1*,*\@igspnih-n3*,*\@igspnih-n4*,*\@igspnih-n6*,*\@igspncbc*
cd /home/ohler/bmajoros/MUMMIE/nicky
hostname
make-42state-site.pl $X.$TRAINING_SEEDS $TRANS_6MER $TRANS_7MER $TRANS_8MER 4 WT-cons.schema 7 bg.hmm $X.site.hmm
echo trace4
cp $SUBMODELS $iter.$sampleIndex.submodels
echo trace5
sub $iter.$sampleIndex.submodels site.hmm $X.site.hmm
echo trace6
MUMMIE/model-combiner $METAMODEL $iter.$sampleIndex.submodels $X.hmm
MUMMIE/hmm-edit $X.hmm DTRK phastcons
MUMMIE/hmm-edit $X.hmm VAR all 0 $VARIANCE
echo trace7
rm $iter.$sampleIndex.submodels
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
$EVAL_SCRIPT $X.predictions.gff WT-groups.gff > $X.mummie.eval
echo trace11
";
    close(OUT);
    system("qsub $X.q >& /dev/null");
  }
}


print "\nall jobs have been submitted\n";

