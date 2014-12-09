#!/usr/bin/perl
use strict;
use GffReader;
use ProgramName;
use FastaReader;

my $name=ProgramName::get();
die "$name <in.gff> <seeds.fasta>\n" unless @ARGV==2;
my ($gffFile,$seedsFile)=@ARGV;

my %motifs;
my $reader=new FastaReader($seedsFile);
while(1) {
  my ($defline,$sequence)=$reader->nextSequence();
  last unless $defline;
  $defline=~/>(\S+)/ || die;
  $motifs{$sequence}->{$1}=1;
}

my $reader=new GffReader();
my $sites=$reader->loadGFF($gffFile);

my $n=@$sites;
for(my $i=0 ; $i<$n ; ++$i) {
  my $site=$sites->[$i];
  my $seq;
  my $extra=$site->{additionalFields};
  foreach my $field (@$extra) { if($field=~/seq=(\S+);/) { $seq=$1 } }
  my $miRNAs=$motifs{$seq};
  my @miRNAs=keys %$miRNAs;
  my @extra=@$extra;
  foreach my $miRNA (@miRNAs) {
    my @ex=@extra;
    push @ex,"miRNA=$miRNA;";
    $site->{additionalFields}=\@ex;
    my $gff=$site->toGff();
    print "$gff";
  }
}


