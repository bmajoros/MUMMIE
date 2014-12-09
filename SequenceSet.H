/****************************************************************
 SequenceSet.H
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_SequenceSet_H
#define INCL_SequenceSet_H
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/Vector.H"
#include "BOOM/LinearFunc.H"
#include "EmissionSequence.H"
using namespace std;
using namespace BOOM;

class SequenceSet {
public:
  SequenceSet();
  void load(const String &directoryPath,Schema &,int maxFiles=-1);
  void save();
  int size() const;
  void resize(int);
  EmissionSequence *operator[](int);
  const EmissionSequence *operator[](int) const;
  void normalize();
  void normalize(const BOOM::Vector<LinearFunc> &);
  void saveTransforms(ofstream &) const;
  void loadTransforms(ifstream &);
  int getDim() const;
  BOOM::Vector<EmissionSequence*> &peek() {return seqs;}
  void recode(int N); // recode to Nmers
  void recode(int trackNum,int N); // recode to Nmers
  void addNoise(double mean,double var);
protected:
  BOOM::Vector<EmissionSequence*> seqs;
  BOOM::Vector<LinearFunc> transforms;
};

#endif

