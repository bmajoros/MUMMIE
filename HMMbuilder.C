/**************************************************************
 HMMbuilder.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/

#include "HMMbuilder.H"
#include "BOOM/HigherOrderAlphabet.H"
#include <gsl/gsl_linalg.h>

HMM *HMMbuilder::randomHMM(int numStates,float transitionDensity,
			   int numComp,const Schema &schema,int order,
			   TransitionList *transList,bool uniformChains)
{
  // Instantiate a new HMM
  HMM &hmm=*new HMM(numStates);
  hmm.getSchema()=schema;
  hmm.changeOrder(order);

  // Set transition probabilities
  if(transList) {
    TransitionList::iterator cur=transList->begin(), end=transList->end();
    for(; cur!=end ; ++cur) {
      StatePair p=*cur;
      const int i=p.from, j=p.to;
      if(j==0) hmm.setTransitionProb(i,j,RandomFloat(0.00001,0.001));
      else if(i==j) hmm.setTransitionProb(i,j,RandomFloat(0.8,0.99));
      else hmm.setTransitionProb(i,j,RandomFloat(0.01,0.2));
    }
  }
  else {
    for(int i=0 ; i<numStates-1 ; ++i) hmm.addTransition(i,i+1);
    hmm.addTransition(numStates-1,0);
    hmm.setTransitionProb(numStates-1,0,RandomFloat(0.00001,0.001));
    for(int i=0 ; i<numStates ; ++i)
      for(int j=0 ; j<numStates ; ++j)
	if(Random0to1()<transitionDensity) {
	  if(j==0) hmm.setTransitionProb(i,j,RandomFloat(0.00001,0.001));
	  else if(i==j) hmm.setTransitionProb(i,j,RandomFloat(0.8,0.99));
	  else hmm.setTransitionProb(i,j,RandomFloat(0.01,0.2));
	}
  }
  hmm.setTransitionProb(0,0,0.0);
  hmm.normalizeTransitions();

  // Set continuous emission distributions
  int numTracks=schema.getNumContinuous();
  for(int q=0 ; q<numStates ; ++q) {
    GaussianMixture &mix=hmm.getEmissionDistr(q);
    mix.resize(numComp);
    for(int i=0 ; i<numComp ; ++i) mix.setCoef(i,Random0to1());
    mix.normalizeCoefs();
  }
  for(int i=0 ; i<numComp ; ++i) {
    GSL::Vector means(numTracks);
    GSL::Matrix R(numTracks,numTracks), cov(numTracks,numTracks);
    for(int j=0 ; j<numTracks ; ++j) {
      means[j]=Random0to1();
      for(int k=0 ; k<numTracks ; ++k) R(j,k)=Random0to1();
    }
    
    // Make cov be positive definite symmetric
    GSL::Vector work(numTracks), S(numTracks);
    GSL::Matrix V(numTracks,numTracks);
    gsl_linalg_SV_decomp(R.peek(),V.peek(),S.peek(),work.peek());
    GSL::Matrix diag(numTracks,numTracks);
    diag.setAllTo(0.0);
    for(int i=0 ; i<numTracks ; ++i) diag(i,i)=Random0to1();
    GSL::Matrix tmp;
    R.times(diag,tmp);
    R.transpose();
    tmp.times(R,cov);
    MultiGauss G(means,cov);
    for(int q=0 ; q<numStates ; ++q) hmm.getEmissionDistr(q).setDistr(i,G);
  }

  // Set discrete emission distributions
  Array2D< Array1D<double> > &discrete=hmm.getDiscreteEmitDistr();
  discrete.resize(numStates,numTracks);
  numTracks=schema.getNumDiscrete();
  for(int q=0 ; q<numStates ; ++q) {
    for(int i=0 ; i<numTracks ; ++i) {
      Array1D<double> &discreteSlice=discrete[q][i];
      Alphabet &alphabet=schema.getAlphabet(i);
      const double uniformP=1.0/alphabet.size();
      HigherOrderAlphabet hoa(alphabet,order+1);
      const int numNmers=hoa.getNumNmers();
      discreteSlice.resize(numNmers);
      for(NmerSymbol nmer=0 ; nmer<numNmers ; ++nmer)
	discreteSlice[nmer]=uniformChains ? uniformP : Random0to1();
    }
  }
  hmm.normalizeDiscreteEmit();

  // Return it
  return &hmm;
}

