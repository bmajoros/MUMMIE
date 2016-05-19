/****************************************************************
 Schema.H
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_Schema_H
#define INCL_Schema_H
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/Vector.H"
#include "BOOM/Alphabet.H"
#include "BOOM/Map.H"
using namespace std;
using namespace BOOM;

class Schema {
public:
  Schema(int numDiscrete=0,int numContinuous=0);
  Schema(const String &filename);
  int getNumDiscrete() const;
  int getNumContinuous() const;
  void setNumDiscrete(int);
  void setNumContinuous(int);
  void load(istream &);
  void save(ostream &) const;
  const String &getDiscreteName(int);
  const String &getContinuousName(int);
  void setDiscreteName(int,const String &);
  void setContinuousName(int,const String &);
  void addContinuousTrack(const String &name);
  void dropContinuousTrack(const String &name);
  void dropDiscreteTrack(const String &name);
  Alphabet &getAlphabet(int);
  int lookupDiscreteID(const String &name); // -1 if not found
  int lookupContinuousID(const String &name);// -1 if not found
protected:
  int numDiscrete, numContinuous;
  Vector<String> discreteNames, continuousNames;
  Map<String,int> discreteIDs, continuousIDs;
  Vector<Alphabet> alphabets;
};

#endif

