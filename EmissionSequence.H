/****************************************************************
 EmissionSequence.H
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_EmissionSequence_H
#define INCL_EmissionSequence_H
#include <iostream>
#include "BOOM/Vector.H"
#include "BOOM/Alphabet.H"
#include "BOOM/Sequence.H"
#include "Emission.H"
#include "Schema.H"
using namespace std;
using namespace BOOM;

class EmissionSequence {
public:
  EmissionSequence(const Schema &);
  void append(const Emission &);
  void append(const EmissionSequence &);
  void recode(int N); // converts discrete tracks to Ngrams
  void recode(int trackNum,int N); // convert to Ngrams
  void unencode(int N); // convert back to 0th order from Ngrams
  int length() const;
  const Emission &operator[](int) const;
  Emission &operator[](int i) {return S[i];}
  int getNumTracks() const;
  Schema &getSchema();
  void save(const String &filename) const;
  void save(ostream &) const;
  void scaleAllTracks(double multiplier);
  void setFilename(const String &);
  const String &getFilename() const;
  void copySeqAndSchema(const EmissionSequence &);
  EmissionSequence *getSubsequence(int begin,int length);
  Sequence *getDiscreteSeq(int trackID);
  void getExtrema(int trackID,float &minVal,float &maxVal);
  void dropDiscreteTrack(int trackID);
protected:
  Vector<Emission> S;
  Schema schema;
  String filename;
};

#endif

