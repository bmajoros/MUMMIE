/**************************************************************
 FastViterbi.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#include <math.h>
#include <float.h>
#include <iostream>
#include "FastViterbi.H"
#include "BOOM/Constants.H"
#include "ForwardAlgorithm.H"
#include "BackwardAlgorithm.H"
using namespace std;

#define MOST_NEGATIVE -DBL_MAX

FastViterbi::FastViterbi(HMM &hmm,bool posterior)
  : numStates(hmm.countStates()), hmmGraph(hmm), posterior(posterior)
{
  // ctor
}



BOOM::Vector<STATE> *FastViterbi::buildPath(int L,STATE finalState,
					    Array2D<short> &linkBack)
{
  int currentState=finalState;
  BOOM::Vector<STATE> &path=*new BOOM::Vector<STATE>(L);
  path[L-1]=finalState;
  for(int pos=L-1 ; pos>0 ; --pos) {
    currentState=linkBack[currentState][pos+1];
    path[pos-1]=currentState;
  }
  return &path;
}



BOOM::Vector<STATE> *FastViterbi::buildPath(int L,STATE finalState,
					    Array2D<short> &linkBack,
					    BOOM::Vector<double> &scores,
					    ForwardAlgorithm &fw,
					    BackwardAlgorithm &bw)
{
  const double likelihood=fw.getLogP();
  int currentState=finalState;
  BOOM::Vector<STATE> &path=*new BOOM::Vector<STATE>(L);
  scores.resize(L);
  path[L-1]=finalState;
  scores[L-1]=fw(finalState,L)+bw(finalState,L)-likelihood;
  for(int pos=L-1 ; pos>0 ; --pos) {
    currentState=linkBack[currentState][pos+1];
    path[pos-1]=currentState;
    //double score=fw(currentState,pos+1)+bw(currentState,pos+1)-likelihood;
    double score=fw(currentState,pos)+bw(currentState,pos)-likelihood;//###
    scores[pos-1]=score;
 }
  return &path;
}



/* This method returns the state sequence of the most probable path (not 
   including the start state). */
BOOM::Vector<STATE> *FastViterbi::getPath(const EmissionSequence &sequence,
					  double &score,BOOM::Vector<double>
					  *stateScores)
{
  ForwardAlgorithm *fw;
  BackwardAlgorithm *bw;
  double likelihood;
  score=0;
  const HMM &hmm=hmmGraph.getHMM();
  if(posterior) {
    fw=new ForwardAlgorithm(hmmGraph,sequence);
    bw=new BackwardAlgorithm(hmmGraph,sequence);
    likelihood=fw->getLogP();
  }
  int L=sequence.length(), Lplus1=L+1;
  Array2D<short> linkBack(numStates,Lplus1);
  Array2D<double> m(numStates,Lplus1);
  m.setAllTo(NEGATIVE_INFINITY);
  m[0][0]=0.0;
  for(int i=1 ; i<=L ; ++i) {
    const Emission &s=sequence[i-1];
    bool foundPred=false;
    for(int q=1 ; q<numStates ; ++q) {
      const BOOM::Vector<StateDoublePair> &precedingStates=
	hmmGraph.statesPreceding(q);
      int numPreceding=precedingStates.size();
      double bestP=NEGATIVE_INFINITY;
      int bestPredecessor;
      for(int j=0 ; j<numPreceding ; ++j) {
	StateDoublePair precedingPair=precedingStates[j];
	STATE y=precedingPair.state;
	double inductiveP=m[y][i-1];
	double newP=posterior ? inductiveP :
	  safeAdd(inductiveP,precedingPair.logP);
	if(newP>bestP){
	  bestP=newP;
	  bestPredecessor=y;
	}
      }
      double emission=posterior ?
	(*fw)(q,i)+(*bw)(q,i)-likelihood :
	hmm.getEmissionProb(q,s);
      m[q][i]=safeAdd(bestP,emission);
      if(isFinite(m[q][i])) foundPred=true;
      linkBack[q][i]=bestPredecessor;
    }
  }
  int bestState=1;
  double bestScore=NEGATIVE_INFINITY;
  for(int i=1 ; i<numStates ; ++i) {
    const double score=posterior ? m[i][L] :
      safeAdd(m[i][L],log(hmm.getTransitionProb(i,0)));
    if(score>bestScore){
      bestState=i;
      bestScore=score;
    }
  }
  if(!isFinite(bestScore)) throw "No path";
  score=bestScore;
  BOOM::Vector<STATE> *retval=posterior && stateScores ? 
    buildPath(L,bestState,linkBack,*stateScores,*fw,*bw) :
    buildPath(L,bestState,linkBack);
  if(posterior) {delete fw; delete bw;}
  return retval;
}



