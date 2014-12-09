/****************************************************************
 BackwardAlgorithm.H
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_BackwardAlgorithm_H
#define INCL_BackwardAlgorithm_H
#include "HMM.H"
#include "BOOM/Sequence.H"
#include "BOOM/Array1D.H"
#include "BOOM/Array2D.H"

/**********************************************************************
 This is "Backward Algorithm."  It computes 
    b(k,i) = P[an HMM M will next emit x(i+1)..x(L) and then terminate
               | M is currently in state k]
 Compare this to the probability computed by the Forward Algorithm.
 The HMM begins in state 0 (which emits nothing) and continues its
 computation until it enters zero again (which again emits nothing).
 *********************************************************************/

class BackwardAlgorithm
{
public:
  BackwardAlgorithm(HMM &,Sequence &);
  double operator()(int state,int sequenceIndex); // Uses 1-based indices
  double getP();
private:
  const int numStates, seqLen;
  Array2D<double> dpMatrix;
  double P, sum;
  HMM &hmm;
  Sequence &sequence;
  void compute();
};


#endif
