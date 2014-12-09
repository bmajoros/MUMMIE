/****************************************************************
 Emission.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "Emission.H"
using namespace std;
using namespace BOOM;


Emission::Emission() 
  : hashCache(UINT_MAX)
{
  // ctor
}



void Emission::resize(int numDiscrete,int numContinuous)
{
  continuous.resize(numContinuous);
  discrete.resize(numDiscrete);
}



int Emission::getDiscreteDim() const
{
  return discrete.size();
}



int Emission::getContinuousDim() const
{
  return continuous.getDim();
}



GSL::Vector &Emission::getContinuous()
{
  return continuous;
}



NmerSymbol &Emission::getDiscrete(int index)
{
  return discrete[index];
}



NmerSymbol Emission::getDiscrete(int index) const
{
  return discrete[index];
}



void Emission::printOn(ostream &os) const
{
  int nd=discrete.size(), nc=continuous.getDim();
  os<<"[";
  for(int i=0 ; i<nd ; ++i) {
    os<<discrete;
    if(i<nd-1 || nc>0) os<<"\t";
  }
  for(int i=0 ; i<nc ; ++i) {
    os<<continuous;
    if(i<nc-1) os<<"\t";
  }
  os<<"]"<<endl;
}



ostream &operator<<(const Emission &e,ostream &os)
{
  e.printOn(os);
  return os;
}



unsigned Emission::hash() const
{
  if(hashCache!=UINT_MAX) return hashCache;
  long unsigned h=0;
  int n=continuous.getDim();
  for(int i=0 ; i<n ; ++i) {
    double c=continuous[i];
    long unsigned uc=*(long unsigned*)(void*)&c;
    h=(h<<1)^uc;
    //h=(h<<8)^uc;
    //h=(unsigned)(int(h)-int(uc));
  }
  //cout<<"h="<<h<<endl;
  n=discrete.size();
  for(int i=0 ; i<n ; ++i) {
    unsigned d=(unsigned)discrete[i];
    h=(h<<1)^d;
  }
  return const_cast<Emission*>(this)->hashCache=(unsigned) h;
}



bool Emission::operator==(const Emission &e) const
{
  if(!(continuous==e.continuous)) return false;
  int n=discrete.size();
  NmerSymbol *s=discrete.getRawArray(), *se=e.discrete.getRawArray();
  for(int i=0 ; i<n ; ++i, ++s, ++se) if(*s!=*se) return false;
  return true;
}





