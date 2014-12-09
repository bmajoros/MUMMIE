/****************************************************************
 BackwardAlgorithm.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <math.h>
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/SumLogProbs.H"
#include "BackwardAlgorithm.H"


BackwardAlgorithm::BackwardAlgorithm(const HMMGraph &graph,
				     const EmissionSequence &sequence)
  : graph(graph),
    hmm(graph.getHMM()),
    dpMatrix(graph.getHMM().countStates(),sequence.length()+1),
    sequence(sequence),
    numStates(graph.getHMM().countStates()),
    seqLen(sequence.length())
{
  compute();
}



double BackwardAlgorithm::getLogP() const
{
  return P;
}



double BackwardAlgorithm::operator()(int s,int i) const
{
  return dpMatrix[s][i];
}



void BackwardAlgorithm::compute()
{
  for(int k=0 ; k<numStates ; ++k) {
    dpMatrix[k][seqLen]=hmm.getLogTransProb(k,0);
    //cout<<"ZZZ "<<k<<" "<<seqLen<<" "<<hmm.getLogTransProb(k,0)<<endl;
  }
  //INTERNAL_ERROR;
  for(int i=seqLen-1 ; i>=0 ; --i)
    for(int k=0 ; k<numStates ; ++k) {
      BOOM::Vector<double> V;
      const BOOM::Vector<StateDoublePair> &next=graph.statesFollowing(k);
      BOOM::Vector<StateDoublePair>::const_iterator cur=next.begin(), 
	end=next.end();
      for(; cur!=end ; ++cur) {
	const StateDoublePair &sdp=*cur;
	STATE l=sdp.state;
	if(l==0) continue;
	double transP=sdp.logP;
	double inductiveP=dpMatrix[l][i+1];
	if(!isFinite(inductiveP)) continue;
	double emitP=hmm.getEmissionProb(l,sequence[i]);
	if(!isFinite(emitP)) continue;
	V.push_back(inductiveP+transP+emitP);
      }
      dpMatrix[k][i]=sumLogProbs(V);
    }
  P=dpMatrix[0][0];
  //cout<<dpMatrix<<endl;
  //INTERNAL_ERROR;
}



void BackwardAlgorithm::printOn(ostream &os) const
{
  os<<dpMatrix<<endl;
}



ostream &operator<<(ostream &os,const BackwardAlgorithm &B)
{
  B.printOn(os);
  return os;
}




