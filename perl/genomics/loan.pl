#!/usr/bin/perl
use strict;

die "loan.pl <amount> <#years> <yearly-interest-rate (integer)>\n" 
    unless @ARGV==3;
my ($P,$years,$yearlyR)=@ARGV;

$yearlyR/=100;
my $n=12*$years;
my $r=$yearlyR/12;

print "P=$P r=$r n=$n\n";
my $payment=int(($P*(($r+1)**$n)*$r)/(($r+1)**$n-1)+0.5);
print "\$$payment/month\n";

