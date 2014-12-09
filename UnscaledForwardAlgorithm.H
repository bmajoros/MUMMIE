/**************************************************************
 UnscaledForwardAlgorithm.H
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#ifndef INCL_ForwardAlgorithm_H
#define INCL_ForwardAlgorithm_H
#include "BOOM/Array2D.H"
#include "BOOM/Array1D.H"
#include "HMM.H"

class ForwardAlgorithm
{
public:
  ForwardAlgorithm(HMM &,Sequence &);
  double getP();
  double operator()(int state,int pos); // uses 1-based indices...
private:
  double P;
  Sequence &sequence;
  int numStates, seqLen;
  Array2D<double> dpMatrix;
  HMM &hmm;
  void computeDPMatrix();
};


#endif
