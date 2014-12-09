#!/usr/bin/perl
use strict;
use LinearRegressor;

my @vars=(1.5, 1.25, 0.75, 0.5, 0.375, 0.25, 0.20, 0.15, 0.1, 0.075,
	  0.05, 0.01, 0.005);
my @SNRs=(9.97, 9.54, 7.1, 6.21, 5.15, 3.95, 3.65, 3.01, 2.62, 2.49,
	  2.32, 2.05, 2.02);
my @sens=(0.05, 0.05, 0.08, 0.14, 0.2, 0.33, 0.4, 0.46, 0.5, 0.51,
	  0.52, 0.53, 0.53);
#my @SNRs=(11.5, 12.19, 13.28, 14.95, 14.43, 10.41, 9.4, 6.52, 4.46,
#	  3.2, 2.75, 2.2, 2.16);
#my @sens=(0.09, 0.09, 0.1, 0.11, 0.12, 0.16, 0.19, 0.25, 0.41, 0.54,
#	  0.58, 0.6, 0.62);
my $N=@vars;
my %types;
for(my $i=3 ; $i<$N ; ++$i) {
  my $var=$vars[$i];
  #my $SNR=$SNRs[$i];
  my $SNR=$sens[$i];
  next unless -e "uniform-var$var-revised.gff";
  my (@counts,$sum);
  open(IN,"cat uniform-var$var-revised.gff | gff-get-extra-field.pl type | sort | uniq -c|");
  while(<IN>) {
    chomp;
    if(/(\d+)\s+(\S+)/) {
      my ($count,$type)=($1,$2);
      $sum+=$count;
      push @counts,[$count,$type];
    }
  }
  close(IN);
  foreach my $pair (@counts) {
    my ($count,$type)=@$pair;
    my $p=$count/$sum;
    push @{$types{$type}},[$SNR,$p];
  }
}
my $regressor=new LinearRegressor();
my @types=keys %types;
foreach my $type (@types) {
  print "$type\t";
  my $array=$types{$type};
  my (@X,@Y);
  foreach my $pair (@$array) {
    my ($snr,$p)=@$pair;
    print "$snr\t$p\n";
    push @X,$snr;
    push @Y,$p;
  }
  print "\n";
  my ($slope,$intercept,$r,$r2)=$regressor->regress(\@X,\@Y);
  print "$r\n";
}
