/**************************************************************
 PosteriorViterbi.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#include <math.h>
#include <float.h>
#include <iostream>
#include "PosteriorViterbi.H"
#include "BOOM/Constants.H"
#include "BOOM/Map.H"
#include "Posterior.H"
using namespace std;

#define MOST_NEGATIVE -DBL_MAX

inline double safeAdd(double a,double b)
{
  if(a==MOST_NEGATIVE || b==MOST_NEGATIVE) return MOST_NEGATIVE;
  else return a+b;
}


PosteriorViterbi::PosteriorViterbi(HMM &hmm,
				   bool shouldAdd)
  : numStates(hmm.countStates()), hmmGraph(hmm),
    addRatherThanMultiply(shouldAdd)
{
  // ctor
}



Vector<int> *PosteriorViterbi::buildPath(int L,int finalState,
					Array2D<short> &linkBack)
{
  int currentState=finalState;
  Vector<int> &path=*new Vector<int>(L);
  path[L-1]=finalState;
  for(int pos=L-1 ; pos>0 ; --pos)
    {
      currentState=linkBack[currentState][pos+1];
      path[pos-1]=currentState;
    }
  return &path;
}



/*
  This method returns the state sequence of
  the most probable path (not including the
  start state).
 */
Vector<int> *PosteriorViterbi::getPath(Sequence &sequence)
{
  /*
  HMM &hmm=hmmGraph.getHMM();
  Set<int> dummySet;
  Posterior posterior(hmmGraph,sequence,dummySet,dummySet,dummySet,0);
  int L=sequence.getLength(), Lplus1=L+1;
  Array2D<short> linkBack(numStates,Lplus1);
  Array2D<double> m(numStates,Lplus1);
  m.setAllTo(NEGATIVE_INFINITY);
  m[0][0]=0.0;
  for(int i=1 ; i<=L ; ++i)
    {
      Symbol s=sequence[i-1];
      Vector<StateDoublePair> &statesEmittingS=
	hmmGraph.statesEmitting(s);
      int numEmitting=statesEmittingS.size();
      for(int l=0 ; l<numEmitting ; ++l)
	{
	  StateDoublePair &currentPair=statesEmittingS[l];
	  Vector<StateDoublePair> &precedingStates=
	    hmmGraph.statesPreceding(currentPair.state);
	  int numPreceding=precedingStates.size();
	  double bestP=log(0.0);
	  int bestPredecessor;
	  for(int j=0 ; j<numPreceding ; ++j)
	    {
	      StateDoublePair precedingPair=precedingStates[j];
	      double inductiveP=m[precedingPair.state][i-1];

	      //#####################################
	      // multiplying by the transition prob....
	      double newP=inductiveP;
		//safeAdd(inductiveP,precedingPair.logP); //################
	      //#####################################

	      if(newP>bestP)
		{
		  bestP=newP;
		  bestPredecessor=precedingPair.state;
		}
	    }
	  //m[currentPair.state][i]=safeAdd(bestP,currentPair.logP);
	  double post=posterior(currentPair.state,i);
	  if(addRatherThanMultiply) post=exp(post);
	  m[currentPair.state][i]=safeAdd(bestP,post);
	  linkBack[currentPair.state][i]=bestPredecessor;
	}
    }

  int bestState=1;
  double bestScore=NEGATIVE_INFINITY;
  for(int i=0 ; i<numStates ; ++i)
    if(m[i][L-1]>bestScore)
      {
	bestState=i;
	bestScore=m[i][L-1];
      }

  return buildPath(L,bestState,linkBack);
  */
}



Vector<int> *PosteriorViterbi::getPath_Unveil(Sequence &sequence,
				     Set<int> &frameshiftStates,
				     Set<int> &exonStates,
				     Set<int> &startCodonStates,
				     int strandDelta,
				     Set<int> &intronStates,
				     Set<int> &intergenicStates)

{
  HMM &hmm=hmmGraph.getHMM();

  //#################################
  Map<int,int> stateToLabel;
  int N=hmm.countStates();
  for(int i=0 ; i<N ; ++i)
    {
      int fw_state=i;
      if(fw_state>strandDelta) fw_state-=strandDelta;
      if(exonStates.isMember(fw_state)) stateToLabel[i]=0;
      else if(intronStates.isMember(fw_state)) stateToLabel[i]=1;
      else stateToLabel[i]=2;
    }
  //#################################


  Posterior posterior(hmmGraph,sequence,exonStates,intronStates,
		      intergenicStates,strandDelta);
  const int L=sequence.getLength();
  Array2D<short> linkBack(numStates,L+1);
  Array2D<double> m(numStates,L+1);
  Array2D<short> frameshifts(numStates,L+1);
  m.setAllTo(NEGATIVE_INFINITY);
  m[0][0]=0.0;
  frameshifts.setAllTo(0); // ### could be made faster
  for(int i=1 ; i<=L ; ++i)
    {
      Symbol s=sequence[i-1];
      Vector<StateDoublePair> &statesEmittingS=
	hmmGraph.statesEmitting(s);
      int numEmitting=statesEmittingS.size();
      for(int l=0 ; l<numEmitting ; ++l)
	{
	  StateDoublePair &currentPair=statesEmittingS[l];
	  Vector<StateDoublePair> &precedingStates=
	    hmmGraph.statesPreceding(currentPair.state);
	  int numPreceding=precedingStates.size();
	  if(!numPreceding) continue;

	  // #### UNVEIL-specific code:
	  int adjustedCurrentState=currentPair.state;
	  if(adjustedCurrentState>=strandDelta)
	    adjustedCurrentState-=strandDelta;
           // ###

	  double bestP=log(0.0);
	  int bestPredecessor=-1;
	  for(int j=0 ; j<numPreceding ; ++j)
	    {
	      StateDoublePair precedingPair=precedingStates[j];

	      // #### UNVEIL-specific code:
	      int adjustedPrecedingState=precedingPair.state;
	      if(adjustedPrecedingState>=strandDelta)
		adjustedPrecedingState-=strandDelta;
	      // ### THIS HACK IS HIGHLY MODEL-DEPENDENT!
	      if(adjustedCurrentState>=30 && adjustedCurrentState<=113 &&
		 (adjustedPrecedingState<30 || adjustedPrecedingState>113) 
		 && frameshifts[precedingPair.state][i-1]%3)
		continue;
	      // ####

	      double inductiveP=m[precedingPair.state][i-1];

	      //========================> TEMPORARY: MULT BY TRANS PROBS
	      double newP=inductiveP;
	      /*
	      double newP;
	      if(stateToLabel[precedingPair.state]==
		 stateToLabel[currentPair.state])
		newP=inductiveP;
	      else
		newP=safeAdd(inductiveP,precedingPair.logP);
	      */
	      //========================================================

	      if(newP>bestP)
		{
		  bestP=newP;
		  bestPredecessor=precedingPair.state;
		}
	    }
	  if(bestPredecessor<0) continue;//### 5/8/03

	  // #### UNVEIL-specific code:
	  frameshifts[currentPair.state][i]=
	    frameshifts[bestPredecessor][i-1];
	  if(frameshiftStates.isMember(adjustedCurrentState))
	    ++frameshifts[currentPair.state][i];
	  // ###

	  double post=posterior(currentPair.state,i);
	  if(addRatherThanMultiply) post=exp(post); // only for OAD
	  
	  /*if(exonStates.isMember(currentPair.state))
	    {
	      post*=0.5; //##############################!!!!!!!!!!!
	      }
	  */


	  m[currentPair.state][i]=safeAdd(bestP,post);

	  //m[currentPair.state][i]=safeAdd(bestP,currentPair.logP);
	  linkBack[currentPair.state][i]=bestPredecessor;
	}
    }
  int bestState=1;
  double bestScore=NEGATIVE_INFINITY;
  for(int i=0 ; i<numStates ; ++i)
    if(m[i][L-1]>bestScore)
      {
	bestState=i;
	bestScore=m[i][L-1];
      }

  return buildPath(L,1 /* WHY NOT bestState? To preclude partials genes.*/,
		   linkBack);
}







