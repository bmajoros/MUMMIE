#!/usr/bin/perl
use strict;
use ProgramName;
use FastaReader;
use TempFilename;

my $name=ProgramName::get();
die "$name <bulge-seqs.fasta> <num-mix-comp> <schema> <order> <out.hmm>\n"
  unless @ARGV==5;
my ($seedsFile,$mix,$schema,$order,$Outfile)=@ARGV;

my $SUBDIR=TempFilename::generate();
system("mkdir $SUBDIR");

# Load seeds
my @seeds;
my $reader=new FastaReader($seedsFile);
while(1) {
  my ($def,$seq)=$reader->nextSequence();
  last unless $def;
  push @seeds,[$def,$seq];
}

# Make canonical model
my $nextSeed=int(rand(1000000));
my $L=7;
my $numStates=$L+1;
my $outfile="$SUBDIR/canonical$L.hmm";
system("bin/random-HMM -c V -u -s $nextSeed $numStates 1 $mix $schema $order $outfile");
die "bin/random-HMM -c V -u -s $nextSeed $numStates 1 $mix $schema $order $outfile" unless -e $outfile;
++$nextSeed;

my $tgfFile="$SUBDIR/tgf.tgf";
system("perl/make-tgf.pl $schema full $tgfFile");

# Train
my $inHMM="$SUBDIR/canonical$L.hmm";
my $outHMM="$SUBDIR/canonical$L-trained.hmm";
my $nmersDir="$SUBDIR/nmers$L";#TempFilename::generate();
my $tieProfile="$SUBDIR/tie$L.profile";#TempFilename::generate();
makeTieProfile($tieProfile,$L);
system("mkdir $nmersDir");
system("perl/fasta-to-fastb.pl $seedsFile $schema $nmersDir");
my $cmd="bin/baum-welch $inHMM $tgfFile $nmersDir 2 $outHMM -R -s $nextSeed -u -t $tieProfile -n 0.00001";
#print "$cmd\n";
system($cmd);
system("rm $tgfFile");
my $pwd=`pwd`;
chomp $pwd;
system("mv $outHMM $Outfile");

# Set mixture parameters
system("bin/hmm-edit $Outfile MIX ALL 0 0 MIX ALL 1 0 MIX ALL 2 0 MIX ALL 3 1");
system("bin/hmm-edit $Outfile MEAN 0 0 1e-7   MEAN 0 1 0.531313");
system("bin/hmm-edit $Outfile MEAN 1 0 0.5    MEAN 1 1 0.531313");
system("bin/hmm-edit $Outfile MEAN 2 0 1.0    MEAN 2 1 0.531313");
system("bin/hmm-edit $Outfile MEAN 3 0 0.5    MEAN 3 1 0.993844");
system("bin/hmm-edit $Outfile COV ALL ALL ALL 0");
system("bin/hmm-edit $Outfile VAR  0 0 0.01   VAR  0 1 1.10313");
system("bin/hmm-edit $Outfile VAR  1 0 0.01   VAR  1 1 1.10313");
system("bin/hmm-edit $Outfile VAR  2 0 0.01   VAR  2 1 1.10313");
system("bin/hmm-edit $Outfile VAR  3 0 0.01   VAR  3 1 1.15067");

# Clean up
system("rm -r $SUBDIR");

#=====================================================================
sub makeTieProfile {
  my ($filename,$L)=@_;
  open(OUT,">$filename") || die $filename;
  print OUT "fix means\nfix covariance_matrix\nfix transitions\n";
  print OUT "fix weights in states ";
  for(my $q=1 ; $q<$L ; ++$q) { print OUT "$q," }
  print OUT "$L\n";
  close(OUT);
}




