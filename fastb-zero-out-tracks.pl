#!/usr/bin/perl
use strict;
use ProgramName;
use Fastb;

my $name=ProgramName::get();
die "$name <in.fastb> <track,name,list> <out.fastb>
    in.fastb and out.fastb can instead be directory names
" unless @ARGV==3;
my ($infile,$trackList,$outfile)=@ARGV;

if($infile=~/\.fastb$/) { zeroOut($infile,$outfile) }
else { # directory
  my $indir=$infile; my $outdir=$outfile;
  my @files=`ls $indir`;
  my $n=@files;
  for(my $i=0 ; $i<$n ; ++$i) {
    my $file=$files[$i]; chomp $file;
    if(/\/(\S+\.fastb$)/) { $file=$1 }
    zeroOut("$indir/$file","$outdir/$file");
  }
}

sub zeroOut
{
  my ($infile,$outfile)=@_;

  print "$infile => $outfile\n"; return;

  my @tracks=split/,/,$trackList;
  my $fastb=new Fastb($infile);
  foreach my $trackName (@tracks) {
    my $track=$fastb->getTrackByName($trackName);
    my $L=$track->getLength();
    my $data=$track->getData();
    for(my $i=0 ; $i<$L ; ++$i) { $data->[$i]=0 }
  }
  $fastb->save($outfile);
}






