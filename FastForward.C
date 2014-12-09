/**************************************************************
 FastForward.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#include "BOOM/Constants.H"
#include <math.h>
#include "BOOM/String.H"
#include "FastForward.H"
#include <iostream>
#include "HMM.H"
using namespace std;


FastForward::FastForward(HMM &hmm,HMMGraph &hmmGraph,
			 Sequence &sequence,
			 int begin,int len)
  : scalingFactors(len+1),sequence(sequence),
    dpMatrix(hmm.countStates(),len+1),
    numStates(hmm.countStates()),hmm(hmm),seqLen(len),
    hmmGraph(hmmGraph),offset(begin)
{
  computeDPMatrix();
}



FastForward::FastForward(HMM &hmm,HMMGraph &hmmGraph,
			 Sequence &sequence)
  : scalingFactors(sequence.getLength()+1),sequence(sequence),
    dpMatrix(hmm.countStates(),sequence.getLength()+1),
    numStates(hmm.countStates()),hmm(hmm),seqLen(sequence.getLength()),
    hmmGraph(hmmGraph), offset(0)
{
  computeDPMatrix();
}



double FastForward::getLogP()
{
  double sum=log(P);
  for(int i=1 ; i<=seqLen ; ++i) sum+=log(scalingFactors[i]);
  return sum;
}



double FastForward::getScaledP()
{
  return P;
}



const DblArray1D &FastForward::getScalingFactors() const
{
  return scalingFactors;
}



double FastForward::operator()(int state,int pos)
{
  return dpMatrix[state][pos];
}



double FastForward::computeScalingFactor(int i)
{
  double scalingFactor=0, sum;
  for(int l=1 ; l<numStates ; ++l)
    {
      sum=0;
      for(int k=0 ; k<numStates ; ++k)
	sum+=hmm.getTransitionProb(k,l)*dpMatrix[k][i-1];
      scalingFactor+=sum*hmm.getEmissionProb(l,sequence[offset+i-1]);
    }
  if(isNaN(scalingFactor) || isInfinity(scalingFactor))
    throw String("scaling value in forward algorithm = ")+scalingFactor;
  /*if(scalingFactor<0 || scalingFactor>1)
    cerr << "scaling value in forward algorithm = "<<scalingFactor<<endl;*/
  return scalingFactor;
}



bool FastForward::PValueExceeds(FastForward &f)
{
  double difference=log(P)-log(f.P);
  for(int i=1 ; i<=seqLen ; ++i)
    {
      double d=log(scalingFactors[i])-log(f.scalingFactors[i]);
      difference+=d;
    }
  if(isInfinity(difference)||isNaN(difference)) 
    cerr << "FastForward: difference=" << difference << endl;
  bool exceeds=(difference>0);
  return exceeds;
}



void FastForward::computeDPMatrix()
{
  P=0;
  double s_sub_i, sum;
  dpMatrix[0][0]=1;
  for(int i=1 ; i<=seqLen ; ++i) dpMatrix[0][i]=0;
  for(int k=1 ; k<numStates ; ++k) dpMatrix[k][0]=0;
  for(int i=1 ; i<=seqLen ; ++i)
    {
      scalingFactors[i]=s_sub_i=computeScalingFactor(i);
      for(int l=1 ; l<numStates ; ++l)
	{
	  Vector<StateDoublePair> &preceding=
	    hmmGraph.statesPreceding(l);
	  sum=0;
	  Vector<StateDoublePair>::iterator cur=preceding.begin(), 
	    end=preceding.end();
	  for(; cur!=end ; ++cur)
	    {
	      int fromState=(*cur).state;
	      double P=(*cur).P;
	      sum+=P*dpMatrix[fromState][i-1];
	    }
	  double v=dpMatrix[l][i]=sum*(hmm.getEmissionProb(l,
				  sequence[offset+i-1])/s_sub_i) ;
	  if(v<0 || v>1 || isNaN(v) || isInfinity(v))
	    throw String("cell (")+l+","+i+")="+v;
	}
    }
  for(int k=1 ; k<numStates ; ++k)
    P+=hmm.getTransitionProb(k,0)*dpMatrix[k][seqLen];
  if(P==0)
    for(int k=1 ; k<numStates ; ++k)
      {
	cout<<"dpMatrix[k][seqLen]="<<dpMatrix[k][seqLen]<<endl;
	cout<<"hmm.getTransitionProb(k,0)="
	    <<hmm.getTransitionProb(k,0)<<endl;
      }
}


