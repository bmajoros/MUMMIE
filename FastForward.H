/**************************************************************
 FastForward.H
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#ifndef INCL_FastForward_H
#define INCL_FastForward_H
#include "HMMGraph.H"
#include "BOOM/Array2D.H"
#include "BOOM/Array1D.H"
#include "HMM.H"
#include "BOOM/Sequence.H"

class HMM;

class FastForward
{
public:
  FastForward(HMM &,HMMGraph &,Sequence &);
  FastForward(HMM &,HMMGraph &,Sequence &,int begin,int len);
  double getLogP();
  const DblArray1D &getScalingFactors() const;
  bool PValueExceeds(FastForward &);
  double operator()(int state,int pos); // uses 1-based indices...
  double getScaledP();
private:
  double P;
  Sequence &sequence;
  int numStates, seqLen;
  DblArray1D scalingFactors;
  Array2D<double> dpMatrix;
  HMM &hmm;
  HMMGraph &hmmGraph;
  const int offset;
  void computeDPMatrix();
  double computeScalingFactor(int);
};


#endif
