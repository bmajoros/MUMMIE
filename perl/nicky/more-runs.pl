#!/usr/bin/perl
use strict;

go(0.25);
go(0.5);
go(0.05);
go(0.01);

sub go {
  my ($var)=@_;
  system("MUMMIE/hmm-edit 100seeds-cons.hmm VAR all 0 $var");
  system("MUMMIE/parse 100seeds-cons.hmm -d fastb-WT-BLS -I 5-45 -p > 100seeds-cons-var$var-PD.gff");
  system("identify-miRNAs.pl 100seeds-cons-var$var-PD.gff 15percent/mirbase/top100.fasta > tmp.1");
  system("MUMMIE/combine-miRNA-predictions.pl tmp.1 > 100seeds-cons-var$var-PD.gff");
  system("chmod a+r 100seeds-cons-var$var-PD.gff");
}

