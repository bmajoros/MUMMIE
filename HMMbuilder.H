/**************************************************************
 HMMbuilder.H
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#ifndef INCL_HMMbuilder_H
#define INCL_HMMbuilder_H
#include "HMM.H"
#include <iostream>
#include "BOOM/Random.H"
using namespace std;

struct StatePair {
  STATE from, to;
  StatePair(STATE f,STATE t) : from(f), to(t) {}
};

typedef BOOM::Vector<StatePair> TransitionList;

class HMMbuilder {
public:
  HMM *randomHMM(int numStates,float transitionDensity,
		 int numMixtureComponents,const Schema &,int order,
		 TransitionList *onlyTheseTransitions=NULL,
		 bool wantUniformChains=false);
};


#endif
