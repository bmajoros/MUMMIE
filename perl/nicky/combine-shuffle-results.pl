#!/usr/bin/perl
use strict;

my $NUM_ITERATIONS=10; # 100
my $EVAL_SCRIPT="eval-SnSp2.pl"; # "eval-SnSp.pl";

for(my $i=0 ; $i<=$NUM_ITERATIONS ; ++$i) {
  system("cat parallel/$i.*.predictions.gff > parallel/$i.predictions.gff");
  system("$EVAL_SCRIPT parallel/$i.predictions.gff WT-groups.gff > parallel/$i.mummie.eval");
}



