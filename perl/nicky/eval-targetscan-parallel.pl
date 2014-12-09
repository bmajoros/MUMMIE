#!/usr/bin/perl
use strict;
use GffReader;
use SummaryStats;
use ProgramName;
$|=1;

my $name=ProgramName::get();
die "$name <seeds-dir-in> <seeds-dir-out> <conserved:Y/N> <NOCONS:0/1>\n" unless @ARGV==4;
my ($inSeeds,$outSeeds,$consYN,$WANT_NON_CONS)=@ARGV;

system("rm -r $outSeeds") if -e $outSeeds;
system("mkdir $outSeeds");
for(my $sh=0 ; $sh<100 ; ++$sh) {
  system("cat $inSeeds/$sh.training-seeds.fasta | grep '>' > $outSeeds/$sh.seeds");
}

MakeDir("targetscan/temp");
MakeDir("targetscan/q");
MakeDir("targetscan/mpirun");
MakeDir("targetscan/out");

my $reader=new GffReader();
for(my $threshold=-2 ; $threshold<=0.6 ; $threshold+=0.1) {
  for(my $sh=0 ; $sh<100 ; ++$sh) {
    my $id="${sh}_$threshold";
    my $qFile="targetscan/q/$id.q";
    open(OUT,">$qFile") || die;
    print OUT "#!/bin/tcsh
#
#\$ -S /bin/tcsh -cwd
#\$ -o /home/ohler/bmajoros/MUMMIE/nicky/targetscan/mpirun/$id.mpirun -j y
#\$ -l highprio
#\$ -N TGscan$sh
#   #\$ -q *\@ohler*,*\@igspnih-n1*,*\@igspnih-n3*,*\@igspnih-n4*,*\@igspnih-n6*,*\@igspncbc*
cd /home/ohler/bmajoros/MUMMIE/nicky
hostname
eval-targetscan-thread.pl $id $outSeeds/$sh.seeds $consYN $WANT_NON_CONS > targetscan/out/$id.out
";
    close(OUT);
    system("qsub $qFile > /dev/null");
  }
  print "All jobs have been submitted.\n";
}



sub MakeDir {
  my ($dir)=@_;
  if(-e $dir) { system("rm -r $dir") }
  system("mkdir -p $dir");
}
