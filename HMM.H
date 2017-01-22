/**************************************************************
 HMM.H
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#ifndef INCL_HMM_H
#define INCL_HMM_H
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/Array2D.H"
#include "BOOM/Array1D.H"
#include "BOOM/Vector.H"
#include "BOOM/GSL/GaussianMixture.H"
#include "BOOM/GSL/TransitionMatrix.H"
#include "BOOM/HashTable.H"
#include "BOOM/Mutex.H"
#include "BOOM/Set.H"
#include "BOOM/HigherOrderAlphabet.H"
#include "Schema.H"
#include "EmissionSequence.H"
#include "State.H"
using namespace std;
using namespace BOOM;
using namespace GSL;

class NegativeInfinityException {};

class HMM {
public:
  HMM(const String &filename);
  HMM(istream &);
  HMM(int numStates);
  HMM(const HMM &);
  virtual ~HMM();
  const HMM &operator=(const HMM &);
  void addStates(int numNewStates);
  Schema &getSchema() {return schema;}
  double getEmissionProb(int inState,const Emission &,
			 bool omitContinuous=false) const;
  double getEmissionProb(int inState,const EmissionSequence &,
			 bool omitContinuous=false) const;
  void setTransitionProb(int fromState,int toState,double p);
  void setEmissionDistr(int inState,const GaussianMixture &);
  void copyEmissionDistr(const HMM &other,int otherState,int thisState);
  GaussianMixture &getEmissionDistr(int q) {return emissionProb[q];}
  Array2D< Array1D<double> > &getDiscreteEmitDistr();
  bool doesTransitionExist(int from,int to) const;
  int countStates() const;
  int getOrder() const {return order;}
  void changeOrder(int o);
  void printOn(ostream &) const;
  double getTransitionProb(int fromState,int toState) const;
  double getLogTransProb(int fromState,int toState) const;
  virtual bool save(const String &filename) const;
  virtual void load(const String &filename);
  virtual bool save(ostream &) const;
  virtual void load(istream &);
  void normalizeTransitions();
  void addTransition(int from,int to);
  double getLogP(EmissionSequence &) const;
  EmissionSequence *sample(BOOM::Vector<STATE> *statePath=NULL,
			   int maxLen=INT_MAX) const;
  EmissionSequence *sampleFromState(STATE,int length) const;
  int numMixtureComponents() const;
  void logifyTransitions();
  void unlogifyTransitions();
  void normalizeJoint(); // debugging
  void normalizeDiscreteEmit();
  void dropContinuousTracks();
  void dropDiscreteTracks();
  void getEquilibrium(GSL::Vector &);
  int getOrder(int state,int trackNum) {return orders[state][trackNum];}
  void addContinuousTrack(const String &name);
  void addDiscreteTrack(const String &name,const Alphabet &,int order);
  void dropContinuousTrack(const String &name);
  void dropDiscreteTrack(const String &name);
  void addMixtureComponent(); // retains parameters of existing components
  const NmerChain &getChain(int trackNum) {return chains[trackNum];}
  Set<int> &getForegroundStates();
  void changeStateOrder(STATE,int order); //all discrete tracks; destructive
protected:
  int numStates;
  int order;
  Schema schema;
  Array2D<int> orders; // [state][trackNum]
  Array2D<double> transitionProb; // [state][state]
  Array1D<GaussianMixture> emissionProb; // [state]
  Array2D< Array1D<double> > discreteEmitProb; // [state][trackNum][nmer]
     // gives conditional prob of last base given the preceding N-1 bases
  Array1D<NmerChain> chains; // [trackNum]
  Set<int> foregroundStates;
  void loadSchema(istream &);
  void saveSchema(ostream &) const;
  void discardInput(istream &,const String &match);
};

ostream &operator<<(ostream &,const HMM &);

#endif
