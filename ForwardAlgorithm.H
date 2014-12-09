/**************************************************************
 ForwardAlgorithm.H
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#ifndef INCL_ForwardAlgorithm_H
#define INCL_ForwardAlgorithm_H
#include "BOOM/Array2D.H"
#include "BOOM/Array1D.H"
#include "EmissionSequence.H"
#include "HMMGraph.H"

class ForwardAlgorithm {
public:
  ForwardAlgorithm(const HMMGraph &,const EmissionSequence &);
  double getLogP() const;
  double operator()(int state,int pos) const; // uses 1-based indices...
  void printOn(ostream &) const;
private:
  double P;
  const EmissionSequence &sequence;
  int numStates, seqLen;
  Array2D<double> dpMatrix;
  const HMM &hmm;
  const HMMGraph &graph;
  void computeDPMatrix();
};

ostream &operator<<(ostream &,const ForwardAlgorithm &);


#endif
