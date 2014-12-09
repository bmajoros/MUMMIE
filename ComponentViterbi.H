/**************************************************************
 ComponentViterbi.H
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#ifndef INCL_ComponentViterbi_H
#define INCL_ComponentViterbi_H

#include <iostream>
#include "BOOM/Set.H"
#include "HMM.H"
#include "BOOM/BitSet.H"
#include "HMMGraph.H"
#include "BOOM/Vector.H"
#include "BOOM/Array2D.H"
using namespace std;


/**************************************************************
 Just like regular Viterbi, except that components of the 
 mixture model are treated separately as mutually exclusive
 events.  This is useful for Viterbi training, where we have
 to ascribe individual output events to specific mixture
 components so we can re-estimate the component weights and
 parameters.
***************************************************************/

typedef int COMPONENT;
typedef pair<STATE,COMPONENT> StateComponentPair;

class ComponentViterbi {
public:
  ComponentViterbi(HMM &);
  BOOM::Vector<StateComponentPair> *getPath(const EmissionSequence &,
					    double &score);
protected:
  int numStates;
  HMMGraph hmmGraph;
  BOOM::Vector<StateComponentPair> *buildPath(int L,int pathEndState,
					      int pathEndComponent,
					      Array2D<short> &ptr,
					      const EmissionSequence &);
  int pickComponent(int state,const Emission &,double &score);
};


#endif
