/**************************************************************
 HMM.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#include <fstream>
#include "UnscaledHMM.H"
#include "BOOM/Random.H"
#include <iostream>
#include "UnscaledForwardAlgorithm.H"
#include "HMMGraph.H"
#include "BOOM/Constants.H"
using namespace std;


HMM::HMM(Alphabet &alphabet,int numStates)
  : numStates(numStates), emissionProb(numStates,alphabet.getNumElements()),
    transitionProb(numStates,numStates), alphabet(alphabet)
{
  emissionProb.setAllTo(0.0);
  transitionProb.setAllTo(0.0);
}



HMM::HMM(const HMM &other)
  : emissionProb(other.emissionProb), numStates(other.numStates), 
    transitionProb(other.transitionProb), alphabet(other.alphabet)
{
}



HMM::HMM(const String &filename,
				     Alphabet &alphabet)
  : emissionProb(0,0), transitionProb(0,0), alphabet(alphabet)
{
  load(filename);
}



HMM::HMM(istream &is,Alphabet &alphabet)
  : transitionProb(0,0), emissionProb(0,0), alphabet(alphabet)
{
  load(is);
}



void HMM::printOn(ostream &os)
{
  os << "Transitions:\n";
  for(int i=0 ; i<numStates ; ++i)
    {
      for(int j=0 ; j<numStates ; ++j)
	os << '\t' << transitionProb[i][j];
      os << '\n';
    }
  int numSymbols=alphabet.getNumElements();
  os << "\nEmissions:\n";
  for(int i=0 ; i<numStates ; ++i)
    {
      for(int j=0 ; j<numSymbols ; ++j)
	os << '\t' << emissionProb[i][j];
      os << '\n';
    }  
  os << endl;
}



ostream &operator<<(ostream &os,HMM &model)
{
  model.printOn(os);
  return os;
}



double HMM::getTransitionProb(int from,int to)
{
  return transitionProb[from][to];
}



void HMM::load(const String &fname)
{
  ifstream is(fname.c_str());
  if(!is.good()) throw String("Can't open ")+fname;
  load(is);
}



void HMM::load(istream &is)
{
  int numSymbols=alphabet.getNumElements();

  // Load number of states and allocate arrays of that size
  is >> numStates;
  emissionProb.resize(numStates,numSymbols);
  transitionProb.resize(numStates,numStates);

  // Read transition  && emission probabilities
  is >> transitionProb >> emissionProb;
}



void HMM::normalizeTransitions()
{
  double sum;
  for(int i=0 ; i<numStates ; ++i)
    {
      sum=0;
      for(int j=0 ; j<numStates ; ++j)
	sum+=transitionProb[i][j];
      if(sum>0)
	for(int j=0 ; j<numStates ; ++j)
	  transitionProb[i][j]/=sum;
    }
}



bool HMM::save(const String &fname)
{
  ofstream os(fname.c_str());
  if(!os.good()) throw String("Can't create ")+fname;
  bool success=save(os);
  return success;
}



bool HMM::save(ostream &os)
{
  os << numStates << '\n' 
     << transitionProb << '\n' 
     << emissionProb << endl;
  return true;
}



void HMM::reverseComp()
{
  double t;
  Symbol C=alphabet.lookup('C');
  Symbol T=alphabet.lookup('T');
  Symbol A=alphabet.lookup('A');
  Symbol G=alphabet.lookup('G');
  swap(G,C);
  swap(A,T);
  normalizeEmissions();
  for(int i=0 ; i<numStates ; ++i)
    for(int j=i+1 ; j<numStates ; ++j)
      {
	t=transitionProb[j][i];
	transitionProb[j][i]=transitionProb[i][j];
	transitionProb[i][j]=t;
      }

  // NOTE:  do NOT call normalizeTransitions() here unless
  //        you really want to ... otherwise, the reverseComp
  //        operation will not be "reversible" (i.e., a
  //        symmetric binary relation).
}



bool HMM::doesTransitionExist(int from,int to)
{
  bool exists=transitionProb[from][to]>0;
  return exists;
}



void HMM::normalizeEmissions()
{
  double sum;
  int numSymbols=alphabet.getNumElements();
  for(int from=0 ; from<numStates ; ++from)
    {
      sum=0;
      for(int s=0 ; s<numSymbols ; ++s)
	sum+=emissionProb[from][s];
      if(sum>0)
	for(int s=0 ; s<numSymbols ; ++s)
	  emissionProb[from][s]/=sum;
    }
}



void HMM::addTransition(int from,int to)
{
  setTransitionProb(from,to,Random0to1());
}



int HMM::countStates()
{
  return numStates;
}



void HMM::swap(Symbol s1,Symbol s2)
{
  double t;
  for(int state=0 ; state<numStates ; ++state)
    {
      t=emissionProb[state][s1];
      emissionProb[state][s1]=emissionProb[state][s2];
      emissionProb[state][s2]=t;
    }
}



void HMM::setTransitionProb(int from,int to,double prob)
{
  transitionProb[from][to]=prob;
}



double HMM::getEmissionProb(int state,Symbol s)
{
  return emissionProb[state][s];
}



Alphabet &HMM::getAlphabet()
{
  return alphabet;
}



double HMM::getP(Sequence &seq)
{
  try
  {
      ForwardAlgorithm f(*this,seq);
      return f.getP();
  }
  catch(NegativeInfinityException)
  {
      return NEGATIVE_INFINITY;
  }
}



void HMM::randomizeEmissionProbs()
{
  double sum;
  int numSymbols=alphabet.getNumElements();
  for(int i=0 ; i<numStates ; ++i)
    {
      sum=0.0;
      for(int j=0 ; j<numSymbols ; ++j)
	{
	  double r=RandomFloat(0.20,0.80);
	  sum+=r;
	  emissionProb[i][j]=r;
	}
      if(sum>0)
	for(int j=0 ; j<numSymbols ; ++j)
	  { emissionProb[i][j]/=sum; }
    }
}



void HMM::setEmissionProb(int state,Symbol s,double prob)
{
  emissionProb[state][s]=prob;
}



