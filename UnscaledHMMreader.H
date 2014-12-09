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
#include "UnscaledHMM.H"
#include "BOOM/Alphabet.H"

/***********************************************************
  0->1 : 0.3
  0->2 : 0.7
  1->1 : 0.9 
  1->3 : 0.1
  ...etc...
  state 1 : A=0.3 T=0.2 C=0.4 g=0.1
  state 2 : ...etc...

  Note that 0 is the start state AND the final state, and
  it does not emit symbols.  The "state" statements are
  optional.  Any state not having one will be initialized
  with random emission probabilities.
***********************************************************/

class HMMreader
{
  Alphabet &alphabet;
  struct Transition
  {
    int from, to;
    float P;
    Transition(int from,int to,float P) : from(from), to(to), P(P) {}
  };
  struct Emission
  {
    int state;
    char symbol;
    float P;
    Emission(int s,char sym,float P) : state(s), symbol(sym), P(P) {}
  };
  Regex transRegex, emitRegex, assignRegex;

public:
  HMMreader(Alphabet &);
  HMM *read(const String &filename);
};


#endif
