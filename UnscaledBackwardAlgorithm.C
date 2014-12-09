/****************************************************************
 UnscaledBackwardAlgorithm.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "UnscaledBackwardAlgorithm.H"
#include <math.h>
#include <iostream>
#include "BOOM/String.H"


BackwardAlgorithm::BackwardAlgorithm(HMM &hmm,Sequence &sequence)
  : hmm(hmm),
    dpMatrix(hmm.countStates(),sequence.getLength()+1),
    sequence(sequence),
    numStates(hmm.countStates()),
    seqLen(sequence.getLength())
{
  compute();
}



double BackwardAlgorithm::getP()
{
  return P;
}



double BackwardAlgorithm::operator()(int s,int i)
{
  return dpMatrix[s][i];
}



void BackwardAlgorithm::compute()
{
  for(int k=0 ; k<numStates ; ++k) 
    dpMatrix[k][seqLen]=hmm.getTransitionProb(k,0);
  for(int i=seqLen-1 ; i>=0 ; --i)
    for(int k=0 ; k<numStates ; ++k)
      {
	sum=0.0;	
	for(int l=1 ; l<numStates ; ++l)
	  sum+=dpMatrix[l][i+1] * hmm.getEmissionProb(l,sequence[i]) 
	    * hmm.getTransitionProb(k,l);
	dpMatrix[k][i]=sum;
        if(isinf(sum)) throw NegativeInfinityException();
      }
  P=dpMatrix[0][0];
}



