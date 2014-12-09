#!/usr/bin/perl
use strict;
use FastaReader;
$|=1;

my $BULGE_TYPE=2;
my $SKIPTRANS=0.0;
my $USE_CLUSTERS=0;
my $VARIANCE=0.5; # 0.04 is the "standard" setting for model N3
my $DASH_G="-I"; # "-g" or "-I";
my $EVAL_SCRIPT="eval-SnSp2.pl"; # "eval-SnSp2.pl";

my $NUM_ITERATIONS=100;
my $SHUFFLED="bulge/type$BULGE_TYPE-group.fasta";
my $UNSHUFFLED="bulge$BULGE_TYPE.fasta";
my $TRAINING_SEEDS="training-seeds.fasta";
my $METAMODEL="bulge-metamodel.txt";
my $SUBMODELS="bulge-submodels.txt";
my $TEST_DIR="fastb-WT-nocons-longest";
my $SCHEMA="WT-cons.schema";
my $STATES="5-11";
my $PEAK_THRESHOLD=0.999;

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
    my ($id,$seq)=@$pair;
    my $array=$shuffled{$id};
    next unless $array;
    my $n=@$array;
    my $index=int(rand($n));
    my $shuff=$array->[$index];
    push @sample,$pair;
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
#\$ -N shf$iter
#\$ -q *\@ohler*,*\@igspnih-n1*,*\@igspnih-n3*,*\@igspnih-n4*,*\@igspnih-n6*,*\@igspncbc*
cd /home/ohler/bmajoros/MUMMIE/nicky
hostname
MUMMIE/motif-scan.pl $X.$TRAINING_SEEDS $TEST_DIR 0 > $X.scan.gff
echo trace1
sites-closest-to-a-peak.pl WT-peaks.gff $X.scan.gff $clusters > $X.proximal.gff
echo trace2
eval-SnSp.pl $X.proximal.gff groups.gff > $X.naive.eval
echo trace3
make-bulge-site.pl $X.$TRAINING_SEEDS 4 WT-cons.schema 6 $X.site.hmm
echo trace4
cp $SUBMODELS $iter.submodels
echo trace5
sub $iter.submodels bulge-site1.hmm $X.site.hmm
echo trace6
MUMMIE/model-combiner $METAMODEL $iter.submodels $X.hmm
MUMMIE/hmm-edit $X.hmm DTRK phastcons
MUMMIE/hmm-edit $X.hmm VAR all 0 $VARIANCE
echo trace7
rm $iter.submodels
echo trace8
MUMMIE/parse $X.hmm $TEST_DIR -d -p $DASH_G $STATES > $X.predictions.gff.tmp
echo trace9
identify-miRNAs.pl $X.predictions.gff.tmp $X.$TRAINING_SEEDS > $X.predictions.gff.tmp2
rm $X.predictions.gff.tmp
#MUMMIE/combine-miRNA-predictions.pl $X.predictions.gff.tmp2 > $X.predictions.gff
cp $X.predictions.gff.tmp2 $X.predictions.gff
rm $X.predictions.gff.tmp2
echo trace10
$EVAL_SCRIPT $X.predictions.gff groups.gff > $X.mummie.eval
echo trace11
";
  close(OUT);

  system("qsub $X.q >& /dev/null");
}


print "\nall jobs have been submitted\n";

