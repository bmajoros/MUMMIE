/**************************************************************
 UnscaledForwardAlgorithm.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#include "BOOM/Constants.H"
#include <math.h>
#include "BOOM/String.H"
#include "UnscaledForwardAlgorithm.H"
#include <iostream>
using namespace std;

class NegativeInfinityException;

ForwardAlgorithm::ForwardAlgorithm(HMM &hmm,
				   Sequence &sequence)
  : sequence(sequence),
    dpMatrix(hmm.countStates(),sequence.getLength()+1),
    numStates(hmm.countStates()),
    hmm(hmm),
    seqLen(sequence.getLength())
{
  computeDPMatrix();
}



double ForwardAlgorithm::getP()
{
    return P;
}



double ForwardAlgorithm::operator()(int state,int pos)
{
  return dpMatrix[state][pos];
}



void ForwardAlgorithm::computeDPMatrix()
{
  P=0;
  double sum;
  dpMatrix[0][0]=1;
  for(int i=1 ; i<=seqLen ; ++i) dpMatrix[0][i]=0;
  for(int k=1 ; k<numStates ; ++k) dpMatrix[k][0]=0;
  for(int i=1 ; i<=seqLen ; ++i)
    {
      for(int l=1 ; l<numStates ; ++l)
	{
	  sum=0;	
	  for(int k=0 ; k<numStates ; ++k) 
	    sum+=hmm.getTransitionProb(k,l)*dpMatrix[k][i-1];
	  double v=dpMatrix[l][i]=sum*(hmm.getEmissionProb(l,sequence[i-1]));
          if(isInfinity(v)) throw NegativeInfinityException();
	}
    }
  for(int k=1 ; k<numStates ; ++k) 
    P+=hmm.getTransitionProb(k,0)*dpMatrix[k][seqLen];
}



