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
#include "BOOM/Sequence.H"
#include "BOOM/Array2D.H"
using namespace std;


class NegativeInfinityException {};


class HMM
{
public:
  HMM(const String &fname,Alphabet &);
  HMM(istream &,Alphabet &);
  HMM(Alphabet &,int numStates);
  HMM(const HMM &);
  double getEmissionProb(int inState,Symbol);
  void setTransitionProb(int fromState,int toState,double p);
  void setEmissionProb(int inState,Symbol,double p);
  bool doesTransitionExist(int from,int to);
  int countStates();
  Alphabet &getAlphabet();
  void printOn(ostream &);
  double getTransitionProb(int fromState,int toState);
  virtual bool save(const String &filename);
  virtual void load(const String &filename);
  virtual bool save(ostream &);
  virtual void load(istream &);
  void reverseComp();
  void normalizeTransitions();
  void normalizeEmissions();
  void addTransition(int from,int to);
  double getP(Sequence &seq);
  void randomizeEmissionProbs();
private:
  Array2D<double> transitionProb;
  int numStates;
  Array2D<double> emissionProb;
  Alphabet &alphabet;

  void swap(Symbol,Symbol);
};

ostream &operator<<(ostream &,HMM &);

#endif
