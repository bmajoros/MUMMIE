/**************************************************************
 FastViterbi.H
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#ifndef INCL_FastViterbi_H
#define INCL_FastViterbi_H
#include <iostream>
#include "BOOM/Set.H"
#include "BOOM/BitSet.H"
#include "BOOM/Vector.H"
#include "BOOM/Array2D.H"
#include "HMM.H"
#include "HMMGraph.H"
#include "ForwardAlgorithm.H"
#include "BackwardAlgorithm.H"
using namespace std;


class FastViterbi {
public:
  FastViterbi(HMM &,bool posterior=false);
  BOOM::Vector<STATE> *getPath(const EmissionSequence &,double &score,
			       BOOM::Vector<double> *stateScores=NULL);
protected:
  bool posterior;
  int numStates;
  HMMGraph hmmGraph;
  BOOM::Vector<STATE> *buildPath(int L,STATE pathEndState,
				 Array2D<short> &ptr);
  BOOM::Vector<STATE> *buildPath(int L,STATE pathEndState,
				 Array2D<short> &ptr,
				 BOOM::Vector<double> &,
				 ForwardAlgorithm &,
				 BackwardAlgorithm &);
};


#endif
