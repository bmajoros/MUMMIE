/**************************************************************
 HMMGraph.H
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#ifndef INCL_HMMGraph_H
#define INCL_HMMGraph_H
#include "BOOM/Vector.H"
#include "HMM.H"
using namespace std;
using namespace BOOM;

struct StateDoublePair {
  short state;
  double logP, P;
  StateDoublePair(int s,double logP,double P) 
    : state(s), logP(logP), P(P) {}
};
bool operator<(const StateDoublePair &p1,const StateDoublePair &p2);

class HMMGraph {
public:
  HMMGraph(const HMM &);
  const BOOM::Vector<StateDoublePair> &statesPreceding(int state) const;
  const BOOM::Vector<StateDoublePair> &statesFollowing(int state) const;
  const HMM &getHMM() const;
  void updateProbs(); // when the HMM changes, call this!
private:
  const HMM &hmm;
  typedef BOOM::Vector< BOOM::Vector<StateDoublePair> > StateLists;
  StateLists predecessors; // state -> predecessors
  StateLists successors;   // state -> successors
  void init();
};


#endif
