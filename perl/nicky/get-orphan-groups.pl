#!/usr/bin/perl
use strict;
use GffReader;
use ProgramName;

my $name=ProgramName::get();
die "$name <sites.gff> <groups.gff>\n" unless @ARGV==2;
my ($sitesFile,$groupsFile)=@ARGV;

my $reader=new GffReader();
my $sitesHash=$reader->hashBySubstrate($sitesFile);
my $groupsHash=$reader->hashBySubstrate($groupsFile);
my @substrates=keys %$groupsHash;
my $numSubstrates=@substrates;
for(my $i=0 ; $i<$numSubstrates ; ++$i) {
  my $substrate=$substrates[$i];
  #print "substrate $substrate\n";
  my $sites=$sitesHash->{$substrate};
  my $groups=$groupsHash->{$substrate};
  foreach my $group (@$groups) {
    #print "CHECKING ";print $group->toGff();
    my $hit=0;
    foreach my $site (@$sites) {
      if($site->overlapsOther($group) && $site->getScore()>=8) {
	$hit=1;
	last
      }
    }
    if(!$hit) { print $group->toGff() }
  }
}


