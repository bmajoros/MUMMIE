#!/usr/bin/perl
use strict;
use SummaryStats;

my @files=`ls targetscan/out`;
my $numFiles=@files;
my %hash;
for(my $i=0 ; $i<$numFiles ; ++$i) {
  my $file=$files[$i];
  chomp $file;
  $file=~/(\S+)_(\S+)/;
  my ($sh,$threshold)=($1,$2);
  my $line=`cat targetscan/out/$file`;
  push @{$hash{$threshold}},$line;
}

my @keys=sort {$a <=> $b} keys %hash;
my $numKeys=@keys;
for(my $i=0 ; $i<$numKeys ; ++$i) {
  my $key=$keys[$i];
  my $lines=$hash{$key};
  my (@SNR,@PPV,@FDR,@Sn,@Sp,@F);
  foreach my $line (@$lines) {
    chomp $line;
    my @fields=split/\s+/,$line;
    my ($SNR,$PPV,$FDR,$Sn,$Sp,$F)=@fields;
    push @SNR,$SNR; push @PPV,$PPV; push @FDR,$FDR;
    push @Sn,$Sn;   push @Sp,$Sp;   push @F,$F;
  }
  next unless @SNR>0 && @Sn>0;

  report(\@SNR); report(\@PPV); report(\@FDR);
  report(\@Sn);  report(\@Sp);  report(\@F);
  print "\n";
}

sub report {
  my ($A)=@_;
  my ($mean,$stddev,$min,$max)=SummaryStats::roundedSummaryStats($A);
  print "$mean\:$stddev\t";
}
