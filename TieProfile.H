/****************************************************************
 TieProfile.H
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_TieProfile_H
#define INCL_TieProfile_H
#include <iostream>
#include "BOOM/Vector.H"
#include "TieEvent.H"
using namespace std;
using namespace BOOM;

/*
Examples of fixing parameters:
  fix weights in states 1,2,3,4,5,6,7,8
  fix weights
  fix means
  fix covariance_matrix
  fix transitions
  fix transitions 1->1 1->2 1->3 2->1 2->2 2->3 3->1 3->2 3->3
  fix chain in state 3

Examples of tying parameters:
  tie chains between states 3,4
  tie means between states 1,3,5,6
  tie covariances between components
  tie transitions 1->2 3->4 5->6
 */

class TieProfile {
public:
  TieProfile();
  TieProfile(const String &filename);
  int numEvents() const;
  const TieEvent &operator[](int) const;
  void load(const String &filename);
protected:
  void parseTie(const String &line,const BOOM::Vector<String> &fields);
  void parseFix(const String &line,const BOOM::Vector<String> &fields);
  void error(const String &line,const String &extra="");
  BOOM::Vector<TieEvent> events;
};

#endif

