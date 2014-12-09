#!/usr/bin/perl
use strict;
use Fastb;
use ProgramName;

my $name=ProgramName::get();
die "$name <fastb-dir>\n" unless @ARGV==1;
my ($dir)=@ARGV;

my @files=`ls $dir`;
my $numFiles=@files;
for(my $i=0 ; $i<$numFiles ; ++$i) {
  my $filename=$files[$i];
  chomp $filename;
  my $fastb=new Fastb("$dir/$filename");
  #my $track=$fastb->getTrackByName("AGO-ReadCount");
  my $track=$fastb->getTrackByName("AGO-Signal");
  my $groups=$track->getNonzeroRegions();
  foreach my $group (@$groups) {
    my ($begin,$end)=@$group;
    my $data=$track->getData();
    my ($max, $sum);
    for(my $pos=$begin ; $pos<$end ; ++$pos) {
      my $x=$data->[$pos];
      if($x>$max) { $max=$x }
      $sum+=$x;
    }
    #print "$max\t$sum\n";
    print "$max\n";
  }

}


