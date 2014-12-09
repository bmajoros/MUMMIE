/****************************************************************
 Emission.H
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_Emission_H
#define INCL_Emission_H
#include <iostream>
#include "BOOM/GSL/Vector.H"
#include "BOOM/Array1D.H"
#include "BOOM/Symbol.H"
using namespace std;
using namespace BOOM;

/**********************************************************************
  Represents a single emission event from a single state.  May contain
  multiple tracks, each of which may be continuous (double) or discrete
  (symbol from some alphabet).
***********************************************************************/
class Emission {
public:
  Emission();
  void resize(int numDiscrete,int numContinuous);
  int getDiscreteDim() const;
  int getContinuousDim() const;
  GSL::Vector &getContinuous();
  double operator[](int index) const {return continuous[index];}
  NmerSymbol &getDiscrete(int index);
  NmerSymbol getDiscrete(int index) const;
  void printOn(ostream &) const;
  unsigned hash() const;
  bool operator==(const Emission &) const;
protected:
  unsigned hashCache;
  GSL::Vector continuous;
  Array1D<NmerSymbol> discrete;
};

ostream &operator<<(const Emission &,ostream &);

#endif

