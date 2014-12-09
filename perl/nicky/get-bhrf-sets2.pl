#!/usr/bin/perl
use strict;
use ProgramName;

my $FASTB_DIR="fastb-WT-cons-longest";#"fastb-WT";

my $POSDIR1="pos1G";
my $POSDIR3="pos3G";
my $NEGDIR="negG";
system("mkdir $POSDIR1") unless -e $POSDIR1;
system("mkdir $POSDIR3") unless -e $POSDIR3;
system("mkdir $NEGDIR") unless -e $NEGDIR;

