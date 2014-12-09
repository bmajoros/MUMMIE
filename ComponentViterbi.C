/**************************************************************
 ComponentViterbi.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#include <math.h>
#include <float.h>
#include <iostream>
#include "ComponentViterbi.H"
#include "BOOM/Constants.H"
using namespace std;

#define MOST_NEGATIVE -DBL_MAX

/*
inline double safeAdd(double a,double b)
{
  if(a==MOST_NEGATIVE || b==MOST_NEGATIVE) return MOST_NEGATIVE;
  else return a+b;
}
*/


ComponentViterbi::ComponentViterbi(HMM &hmm)
  : numStates(hmm.countStates()), hmmGraph(hmm)
{
  // ctor
}



BOOM::Vector<StateComponentPair> *
ComponentViterbi::buildPath(int L,
			    int finalState,
			    int finalComponent,
			    Array2D<short> &linkBack,
			    const EmissionSequence &s)
{
  int currentState=finalState;
  BOOM::Vector<StateComponentPair> &path=
    *new BOOM::Vector<StateComponentPair>(L);
  path[L-1]=StateComponentPair(finalState,finalComponent);
  double dummy;
  for(int pos=L-1 ; pos>0 ; --pos) {
    currentState=linkBack[currentState][pos+1];
    int comp=pickComponent(currentState,s[pos],dummy);
    path[pos-1]=StateComponentPair(currentState,comp);
  }
  return &path;
}



int ComponentViterbi::pickComponent(int state,const Emission &s,
				    double &score)
{
  GaussianMixture &distr=hmmGraph.getHMM().getEmissionDistr(state);
  return distr.mostProbableComponent(s.getContinuous(),score);
}



/* This method returns the state sequence of the most probable path (not 
   including the start state). */
BOOM::Vector<StateComponentPair> *
ComponentViterbi::getPath(const EmissionSequence &sequence,
			  double &score)
{
  score=0;
  HMM &hmm=hmmGraph.getHMM();
  int L=sequence.length(), Lplus1=L+1;
  Array2D<short> linkBack(numStates,Lplus1);
  Array2D<double> m(numStates,Lplus1);
  m.setAllTo(NEGATIVE_INFINITY);
  m[0][0]=0.0;
  for(int i=1 ; i<=L ; ++i) {
    const Emission &s=sequence[i-1];
    bool foundPred=false;
    for(int q=1 ; q<numStates ; ++q) {
      BOOM::Vector<StateDoublePair> &precedingStates=
	hmmGraph.statesPreceding(q);
      int numPreceding=precedingStates.size();
      //      cout<<"XXX "<<q<<" "<<precedingStates[0].logP<<" "<<precedingStates[1].logP<<" "<<precedingStates[2].logP<<endl;
      //if(numPreceding>0) foundPred=true;
      double bestP=log(0.0);
      int bestPredecessor;
      for(int j=0 ; j<numPreceding ; ++j) {
	StateDoublePair precedingPair=precedingStates[j];
	double inductiveP=m[precedingPair.state][i-1];
	//cout<<inductiveP<<endl;
	double newP=safeAdd(inductiveP,precedingPair.logP);
	if(newP>bestP){
	  bestP=newP;
	  bestPredecessor=precedingPair.state;
	}
      }
      double emission;
      pickComponent(q,s,emission);
      //cout<<"XXX "<<q<<" "<<emission<<endl;
      m[q][i]=safeAdd(bestP,log(emission));
      if(isFinite(m[q][i])) foundPred=true;
      linkBack[q][i]=bestPredecessor;
    }
    /*
    if(!foundPred) {
      for(int q=1 ; q<numStates; ++q)
	cout<<q<<" "<<hmm.getEmissionProb(q,s)<<" "<<s<<" "<<*hmm.getEmissionDistr(q)<<" "<<m[q][i-1]<<endl;
      throw String("No linkback for position ")+(i-1)+", value="+s;
    }
    */
  }
  //cout<<m<<endl;
  //cout<<linkBack<<endl;
  int bestState=1;
  double bestScore=NEGATIVE_INFINITY, dummy;
  for(int i=0 ; i<numStates ; ++i)
    if(m[i][L]>bestScore){
      bestState=i;
      bestScore=m[i][L];
    }
  score=bestScore;
  int bestComp=pickComponent(bestState,sequence[L-1],dummy); // ###
  return buildPath(L,bestState,bestComp,linkBack,sequence);
}



