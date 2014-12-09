/****************************************************************
 IntervalDecoder.H
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_IntervalDecoder_H
#define INCL_IntervalDecoder_H
#include <iostream>
#include "BOOM/GffFeature.H"
#include "BOOM/Array1D.H"
#include "BOOM/Vector.H"
#include "HMM.H"
#include "HMMGraph.H"
#include "EmissionSequence.H"
#include "ForwardAlgorithm.H"
#include "BackwardAlgorithm.H"
using namespace std;
using namespace BOOM;


class IntervalDecoder {
  HMM &hmm;
  int numStates;
  HMMGraph hmmGraph;
  ForwardAlgorithm *fw;
  BackwardAlgorithm *bw;
  Set<int> &FG;
  Array1D<bool> isFg;
  Array1D< Set<int> > activeAtPos; // [pos] -> states w/nonzero posterior there
public:
  IntervalDecoder(HMM &);
  BOOM::Vector<GffFeature> *decode(EmissionSequence &,
				   const String &substrateID, 
				   const String &sourceID, 
				   const String &featureType);
  double getLogP(EmissionSequence &,int begin,int end,ForwardAlgorithm *,
		 BackwardAlgorithm *); // half-open zero-based coordinates!
};

#endif

