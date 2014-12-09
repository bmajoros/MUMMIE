/**************************************************************
 HMMGraph.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/

#include "HMMGraph.H"
#include <iostream>
#include <math.h>
using namespace std;

bool operator<(const StateDoublePair &p1,const StateDoublePair &p2)
{
  return p1.state<p2.state;
}


HMMGraph::HMMGraph(const HMM &hmm)
  : hmm(hmm), 
    predecessors(hmm.countStates()),
    successors(hmm.countStates())
{
  init();
}



void HMMGraph::init()
{
  int numStates=hmm.countStates();

  for(int i=0 ; i<numStates ; ++i)
    for(int j=0 ; j<numStates ; ++j)
      {
	double p=hmm.getTransitionProb(i,j);
	if(p>0)
	  {
	    predecessors[j].push_back(StateDoublePair(i,log(p),p));
	    successors[i].push_back(StateDoublePair(j,log(p),p));
	  }
      }
}



const BOOM::Vector<StateDoublePair> &HMMGraph::statesPreceding(int state) const
{
  return predecessors[state];
}



const BOOM::Vector<StateDoublePair> &HMMGraph::statesFollowing(int state) const
{
  return successors[state];
}



const HMM &HMMGraph::getHMM() const
{
  return hmm;
}



void HMMGraph::updateProbs()
{
  const int Nq=predecessors.size();
  for(STATE q=0 ; q<Nq ; ++q) {
    BOOM::Vector<StateDoublePair> &P=predecessors[q], &S=successors[q];
    BOOM::Vector<StateDoublePair>::iterator cur=P.begin(), end=P.end();
    for(; cur!=end ; ++cur) {
      StateDoublePair &pair=*cur;
      pair.P=hmm.getTransitionProb(pair.state,q);
      pair.logP=log(pair.P);
    }
    cur=S.begin(); end=S.end();
    for(; cur!=end ; ++cur) {
      StateDoublePair &pair=*cur;
      pair.P=hmm.getTransitionProb(q,pair.state);
      pair.logP=log(pair.P);
    }
  }
}


