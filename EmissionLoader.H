/****************************************************************
 EmissionLoader.H
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_EmissionLoader_H
#define INCL_EmissionLoader_H
#include <iostream>
#include "BOOM/Vector.H"
#include "BOOM/String.H"
#include "BOOM/Alphabet.H"
#include "BOOM/Sequence.H"
#include "EmissionSequence.H"
#include "Schema.H"
using namespace std;
using namespace BOOM;


class EmissionLoader {
public:
  EmissionLoader(Schema &);
  EmissionSequence *load(const String &fastbFile);
  void loadContinuous(const String &fastcFile);
  void loadDiscrete(const String &fastaFile,Alphabet &);
  void loadMixed(const String &fastbFile);
  EmissionSequence *compile();
  BOOM::Vector<int> *loadPathFile(const String &filename);
protected:
  BOOM::Vector<Sequence> discrete;
  BOOM::Vector< BOOM::Vector<double> > continuous;
  BOOM::Vector<String> discreteNames, continuousNames;
  int length;
  const Schema schema;
};

#endif

