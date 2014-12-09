#!/usr/bin/perl

open(IN,"/home/bmajoros/bin/colors.txt") ||
  die "can't open colors.txt";
while(<IN>)
  {
    chop;
    push @commands,$_;
  }
close(IN);

$n=@commands;
$r=int(rand($n));
$cmd=$commands[$r];
system($cmd);
print $cmd."\n";

