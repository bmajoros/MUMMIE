/****************************************************************
 Posterior.H
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_Posterior_H
#define INCL_Posterior_H
#include "ForwardAlgorithm.H"
#include "BackwardAlgorithm.H"
#include "HMMGraph.H"
#include "BOOM/Set.H"
#include "BOOM/Map.H"

class Posterior
{
public:
  Posterior(HMMGraph &,Sequence &,Set<int> &exonStates,
	    Set<int> &intronStates,Set<int> &intergenicStates,
	    int strandDelta);
  double operator()(int state,int pos); // uses 1-based indices!
                                        // and returns log scores
private:
  HMM &hmm;
  Sequence &S;
  Array2D<double> dpMatrix;
  HMMGraph &graph;
  Set<int> &exonStates, &intronStates, &intergenicStates;
  Map<int,int> stateToLabel;

  void compute();
};


#endif
