#!/usr/bin/perl
use strict;
use ProgramName;

# Process command line
my $name=ProgramName::get();
die "$name <name>\n" unless @ARGV==1;
my ($name)=@ARGV;

my $Cfile="$name.C";
my $Hfile="$name.H";
die "$Cfile already exists\n" if -s $Cfile;

system("mkdir obj") unless -e "obj";
system("ln -s /home/bmajoros/BOOM") unless -e "BOOM";


###################################################################
#                         Write makefile
###################################################################

my @classes;
if(-e "makefile")
  {
    my $filename="makefile";
    open(MAKEFILE,">>$filename") || 
      die "Can't open $filename for appending";

    # Compile the source file
    print MAKEFILE
"
#--------------------------------------------------------
\$(OBJ)/$name.o:\\
\t\t$name.C \\
\t\t$name.H
\t\$(CC) \$(CFLAGS) -o \$(OBJ)/$name.o -c \\
\t\t$name.C
#---------------------------------------------------------
";

    close(MAKEFILE);
  }


###################################################################
#                         Write .C file
###################################################################

open(CFILE,">$Cfile") || die "can't create $Cfile\n";
print CFILE
"/****************************************************************
 $Cfile
 Copyright (C)2012 William H. Majoros (martiandna\@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include \"$name.H\"
using namespace std;
using namespace BOOM;

$name::$name()
{
  // ctor
}


";
close(CFILE);


###################################################################
#                         Write .H file
###################################################################

open(HFILE,">$Hfile") || die "can't create $Hfile\n";
print HFILE
"/****************************************************************
 $Hfile
 Copyright (C)2012 William H. Majoros (martiandna\@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_RateMatrix_H
#define INCL_RateMatrix_H
#include <iostream>
using namespace std;
using namespace BOOM;


class $name
{
public:
  $name();
};


#endif

";
close(HFILE);

