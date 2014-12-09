#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <in.gff>\n" unless @ARGV==1;
my ($infile)=@ARGV;

my @list;
open(IN,$infile) || die "can't open file: $infile\n";
while(<IN>) {
  chomp;
  next if/^#/;
  my @fields=split/\s+/,$_;
  my $extra=$fields[8];
  my @extra=split/;/,$extra;
  my $phast;
  foreach my $pair (@extra) {
    $pair=~/(\S+)=(\S+)/ || die "can't parse extra field: $pair\n";
    my $field=$1;
    my $value=$2;
    if($field=~/phastcons/) {
      $phast=$value;
      last;
    }
  }
  push @list,[$phast,$_];
}
close(IN);

@list=sort {$b->[0] <=> $a->[0]} @list;
my $n=@list;
for(my $i=0 ; $i<$n ; ++$i) {
  my $pair=$list[$i];
  my ($score,$line)=@$pair;
  print "$line\n";
}
