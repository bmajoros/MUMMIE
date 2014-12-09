/**************************************************************
 HMMreader.H
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#ifndef INCL_HMMreader_H
#define INCL_HMMreader_H

#include "BOOM/String.H"
#include "BOOM/Regex.H"
#include "HMM.H"

/***********************************************************
  0->1 : 0.3
  0->2 : 0.7
  1->1 : 0.9 
  1->3 : 0.1
  ...etc...

  Note that 0 is the start state AND the final state, and
  it does not emit symbols. 
***********************************************************/

class HMMreader
{
  struct Transition {
    int from, to;
    float P;
    Transition(int from,int to,float P) : from(from), to(to), P(P) {}
  };
  Regex transRegex;

public:
  HMMreader();
  HMM *read(const String &filename);
};


#endif
