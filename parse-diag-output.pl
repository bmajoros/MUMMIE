#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use SummaryStats;

my @records;
my $dim;
while(<STDIN>) {
  chomp;
  if(/MEAN=\d+\s+(\S.*\S)/) {
    my @fields=split/\s+/,$1;
    $dim=@fields;
    my $rec=[];
    @$rec=@fields;
    push @records,$rec;
  }
  elsif(/LL=(\S+)/) {
    my $LL=$1;
    my @correlations;
    for(my $i=0 ; $i<$dim ; ++$i) {
      for(my $j=$i+1 ; $j<$dim ; ++$j) {
	my (@X,@Y);
	foreach my $record (@records) {
	  push @X,$record->[$i];
	  push @Y,$record->[$j];
	}
	my $r=SummaryStats::correlation(\@X,\@Y);
	push @correlations,$r;
      }
    }
    my ($aveCor,$sd)=SummaryStats::summaryStats(\@correlations);
    undef @records;
    print "$LL $aveCor\n";
  }
  elsif(/performing EM/){
    my @correlations;
    for(my $i=0 ; $i<$dim ; ++$i) {
      for(my $j=$i+1 ; $j<$dim ; ++$j) {
	my (@X,@Y);
	foreach my $record (@records) {
	  push @X,$record->[$i];
	  push @Y,$record->[$j];
	}
	my $r=SummaryStats::correlation(\@X,\@Y);
	push @correlations,$r;
      }
    }
    my ($aveCor,$sd)=SummaryStats::summaryStats(\@correlations);
    undef @records;
    print "-50000 $aveCor\n";
  }
}


