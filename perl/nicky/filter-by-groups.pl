#!/usr/bin/perl
use strict;
use ProgramName;
use GffReader;

my $name=ProgramName::get();
die "cat file.gff | $name <groups.gff> | ...\n" unless @ARGV==1;
my ($groupsFile)=@ARGV;

my $reader=new GffReader();
my $bySubstrate=$reader->hashBySubstrate($groupsFile);

while(<STDIN>) {
  chomp;
  next if(/^\s*#/);
  my $gff=$_;
  my @fields=split/\s+/,$_;
  next unless @fields>=8;
  my ($substrate,$source,$type,$begin,$end,$score,$strand,$frame,$extra)=
    @fields;
  my $groups=$bySubstrate->{$substrate};
  foreach my $group (@$groups) {
    if($begin<$group->getEnd() && $group->getBegin()<$end) {
      print "$gff\n";
    }
  }
}
