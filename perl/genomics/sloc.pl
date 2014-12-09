#!/usr/bin/perl
use strict;

my $code=`cat *.[HC]`;
my @lines=split/\n/,$code;
my $lines=0;
my $inComment=0;
foreach my $line (@lines)
  {
    while($line=~/(.*)\/\//){$line=$1}
    $line=~s/\/\*.*\*\// /g;
    if($line=~/(.*)\/\*/)
      {
	$inComment=1;
	if($1=~/\S/) {++$lines}
	next;
      }
    if($line=~/\*\/(.*)/)
      {
	$inComment=0;
	$line=$1;
      }
    if($line=~/\S/ && !$inComment) {++$lines}
  }
print "$lines lines\n";


