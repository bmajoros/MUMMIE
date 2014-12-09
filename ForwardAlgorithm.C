/**************************************************************
 ForwardAlgorithm.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#include <iostream>
#include <math.h>
#include "BOOM/String.H"
#include "BOOM/SumLogProbs.H"
#include "BOOM/Constants.H"
#include "ForwardAlgorithm.H"
using namespace std;

class NegativeInfinityException;

ForwardAlgorithm::ForwardAlgorithm(const HMMGraph &graph,
				   const EmissionSequence &sequence)
  : graph(graph), hmm(graph.getHMM()), sequence(sequence),
    dpMatrix(graph.getHMM().countStates(), sequence.length()+1),
    numStates(graph.getHMM().countStates()), seqLen(sequence.length())
{
  computeDPMatrix();
}



double ForwardAlgorithm::getLogP() const
{
  return P;
}



double ForwardAlgorithm::operator()(int state,int pos) const
{
  return dpMatrix[state][pos];
}



void ForwardAlgorithm::computeDPMatrix()
{
  P=0;
  dpMatrix[0][0]=0.0;
  for(int i=1 ; i<=seqLen ; ++i) dpMatrix[0][i]=NEGATIVE_INFINITY;
  for(int k=1 ; k<numStates ; ++k) dpMatrix[k][0]=NEGATIVE_INFINITY;
  for(int i=1 ; i<=seqLen ; ++i) {
    for(int l=1 ; l<numStates ; ++l) {
      BOOM::Vector<double> V;
      const BOOM::Vector<StateDoublePair> &pred=graph.statesPreceding(l);
      BOOM::Vector<StateDoublePair>::const_iterator cur=pred.begin(), 
	end=pred.end();
      for(; cur!=end ; ++cur) {
	const StateDoublePair &sdp=*cur;
	STATE k=sdp.state;
	double transP=sdp.logP;
	double inductiveP=dpMatrix[k][i-1];
	if(!isFinite(inductiveP)) continue;
	V.push_back(inductiveP+transP);
	//cout<<"i="<<i<<" q1="<<l<<" q2="<<k<<" trans="<<transP<<" induct="<<inductiveP<<endl;
      }
      double sum=sumLogProbs(V);
      double emitP=hmm.getEmissionProb(l,sequence[i-1]);
      //cout<<"sum="<<sum<<" emitP="<<emitP<<endl;
      double v=dpMatrix[l][i]=sum+emitP;
      if(isNaN(v)) {
	cout<<sum<<" "<<hmm.getEmissionProb(l,sequence[i-1])<<" "<<endl;
	if(v<0) cout<<v<<"<0"<<endl;
	if(v>1) cout<<v<<">1"<<endl;
	if(isNaN(v)) cout<<v<<" isNaN"<<endl;
	if(isInfinity(v)) cout<<v<<" isInf"<<endl;
	throw String("+cell (")+l+","+i+")="+v;
      }
    }
  }
  BOOM::Vector<double> V;

  const BOOM::Vector<StateDoublePair> &pred=graph.statesPreceding(0);
  BOOM::Vector<StateDoublePair>::const_iterator cur=pred.begin(), 
    end=pred.end();
  for(; cur!=end ; ++cur) {
    const StateDoublePair &sdp=*cur;
    STATE k=sdp.state;
    if(k==0) continue;
    double transP=sdp.logP;
    double inductiveP=dpMatrix[k][seqLen];
    if(!isFinite(inductiveP)) continue;
    V.push_back(inductiveP+transP);
  }
  P=sumLogProbs(V);
  if(!isFinite(P))
    for(int k=1 ; k<numStates ; ++k) {
      cout<<"dpMatrix[k][seqLen]="<<dpMatrix[k][seqLen]<<endl;
      cout<<"hmm.getTransitionProb(k,0)="
	  <<hmm.getTransitionProb(k,0)<<endl;
    }
}



void ForwardAlgorithm::printOn(ostream &os) const
{
  os<<dpMatrix<<endl;
}



ostream &operator<<(ostream &os,const ForwardAlgorithm &B)
{
  B.printOn(os);
  return os;
}



