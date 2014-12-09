/**************************************************************
 PosteriorViterbi.H
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#ifndef INCL_PosteriorViterbi_H
#define INCL_PosteriorViterbi_H

#include <iostream>
#include "BOOM/Set.H"
#include "HMM.H"
#include "BOOM/Sequence.H"
#include "BOOM/BitSet.H"
#include "HMMGraph.H"
#include "BOOM/Vector.H"
#include "BOOM/Array2D.H"
using namespace std;


class PosteriorViterbi
{
public:
  PosteriorViterbi(HMM &,bool addRatherThanMultiply);
  Vector<int> *getPath(Sequence &);
  Vector<int> *getPath_Unveil(Sequence &,
				  Set<int> &frameshiftStates,
				  Set<int> &exonStates,
				  Set<int> &startCodonStates,
				  int strandDelta,
				  Set<int> &intronStates,
				  Set<int> &intergenicStates);
private:
  int numStates;
  HMMGraph hmmGraph;
  bool addRatherThanMultiply;
  Vector<int> *buildPath(int L,int pathEndState,
			     Array2D<short> &ptr);
};


#endif
