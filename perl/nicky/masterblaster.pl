#!/usr/bin/perl
use strict;
$|=1;

system("date");
print "Evaluating MUMMIE...\n";
system("run");
waitForJobs("shufex");

exit; ###

print "Evaluating targetscan...\n";
system("eval-targetscan-parallel.pl fc/shuffling/parallel-var1.25 seeds Y 0");
waitForJobs("TGscan");
system("eval-targetscan-finish.pl > tg.eval.sanity.cons");
system("eval-targetscan-parallel.pl fc/shuffling/parallel-var1.25 seeds N 0");
waitForJobs("TGscan");
system("eval-targetscan-finish.pl > tg.eval.sanity.nocons");
system("eval-targetscan-parallel.pl fc/shuffling/parallel-var1.25 seeds N 1");
waitForJobs("TGscan");
system("eval-targetscan-finish.pl > tg.eval.sanity.NONcons");
system("date");

sub waitForJobs {
  my ($label)=@_;
  while(1) {
    sleep(30);
    open(IN,"qstat -u bmajoros|") || die;
    <IN>;<IN>;
    my $found=0;
    while(<IN>) {
      chomp;
      $_=~/(\S.*\S)/;
      $_=$1;
      my @fields=split/\s+/,$_;
      my ($jobID,$prior,$name,$user,$state,$date,$time,$queue)=@fields;
      next if($name eq "QRLOGIN");
      next if($name eq "masterblas");
      next if($state eq "dt" || $state eq "dr");
      if($name=~/$label/) { ++$found }
    }
    close(IN);
    my $date=`date`;
    chomp $date;
    print "$date\t$found $label jobs remaining\n";
    last unless $found;
  }
}




