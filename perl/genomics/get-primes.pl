#!/usr/bin/perl

die "$0 <max-value>\n" unless @ARGV==1;
($maxValue)=@ARGV;

for($i=2 ; $i<=$maxValue ; ++$i)
{
    for($j=2 ; $i*$j<=$maxValue ; ++$j)
    {
	$product=$i*$j;
	$divisible{$product}=1;
    }
}

for($i=2 ; $i<=$maxValue ; ++$i)
{
    push(@primes,$i) unless $divisible{$i};
}

@primes=sort {$b-$a} @primes;
$n=(@primes<5 ? 0+@primes : 5);
for($i=0 ; $i<$n ; ++$i)
{
    print "$primes[$i]\n";
}


