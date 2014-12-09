#!/usr/bin/perl
use strict;

my @vars=(1.5, 1.25, 0.75, 0.5, 0.375, 0.25, 0.20, 0.15, 0.1, 0.075,
	  0.05, 0.01, 0.005);

foreach my $var (@vars) {
  system("MUMMIE/hmm-edit uniform.hmm VAR all 0 $var");
  system("MUMMIE/parse uniform.hmm -d WT-full -g 5-45 -p > uniform-var$var.paths");
  system("identify-miRNAs.pl uniform-var$var.paths 15percent/mirbase/top50.fasta > uniform-var$var.gff");
  system("MUMMIE/combine-miRNA-predictions.pl uniform-var$var.gff > uniform-var$var-uniq.gff");

  system("cat uniform-var$var.gff | revise-1A.pl > uniform-var$var-revised.gff");

next;

  my ($total,$Y);
  open(IN,"cat uniform-var$var.gff | revise-1A.pl | sort | uniq -c|") ||die;
  while(<IN>) {
    if(/(\d+)\s+Y/){$Y=$1}
  }
  close(IN);
  open(IN,"cat uniform-var$var.gff | gff-get-extra-field.pl type | sort | uniq -c | grep 8mer-m1|") ||die;
  while(<IN>) {
    if(/(\d+)\s+8mer-m1/) {
      $total=$1;
    }
  }
  close(IN);
  $total+=$Y;
  my $percent=$Y/$total;
  print "$percent\n";
}



