#!/usr/bin/perl
use strict;
use SummaryStats;

process("SNR");
process("Sp");
process("Sn");
process("F");
process("PPV");
process("FDR");
process("TP",2000);
process("FP",2000);
system("xgraph TP.txt FP.txt -lw 3 &");
system("xgraph SNR.txt F.txt Sn.txt Sp.txt PPV.txt FDR.txt -lw 3 &");



sub process {
  my ($feature,$divisor)=@_;
  if(!$divisor) {$divisor=1}
  my %hash;
  open(IN,"cat parallel/*.mummie.eval |") || die;
  while(<IN>) {
    if(/threshold=(\S+).*$feature\=(\S+)/) { push @{$hash{$1}},$2/$divisor }
  }
  close(IN);
  my @thresholds=keys %hash;
  my $n=@thresholds;
  #print "$n\n";
  my $outfile="$feature.txt";
  open(OUT,">$outfile") || die;
  @thresholds=sort {$a<=>$b} @thresholds;
  foreach my $threshold (@thresholds) {
    my ($mean,$stddev,$min,$max)=
      SummaryStats::summaryStats($hash{$threshold});
    print OUT "$threshold\t$mean\n";
  }
  close(OUT);
}

