/****************************************************************
 EmissionSequence.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/Constants.H"
#include "EmissionSequence.H"
#include "BOOM/HigherOrderAlphabet.H"
using namespace std;
using namespace BOOM;

EmissionSequence::EmissionSequence(const Schema &schema)
  : schema(schema)
{
  // ctor
}



Schema &EmissionSequence::getSchema()
{
  return schema;
}



void EmissionSequence::append(const Emission &e)
{
  S.push_back(e);
}



int EmissionSequence::length() const
{
  return S.size();
}



const Emission &EmissionSequence::operator[](int i) const
{
  return S[i];
}



int EmissionSequence::getNumTracks() const
{
  return schema.getNumDiscrete()+schema.getNumContinuous();
}



void EmissionSequence::save(const String &filename) const
{
  ofstream os(filename.c_str());
  save(os);
}



void EmissionSequence::save(ostream &os) const
{
  // Write discrete tracks
  int n=schema.getNumDiscrete();
  for(int i=0 ; i<n ; ++i) {
    const Alphabet &alpha=schema.getAlphabet(i);
    const String &name=schema.getDiscreteName(i);
    os<<">"<<name<<endl;
    int col=0;
    Vector<Emission>::const_iterator cur=S.begin(), end=S.end();
    for(; cur!=end ; ++cur) {
      const Emission e=*cur;
      Symbol s=e.getDiscrete(i);
      //os<<int(s)<<' ';
      os<<alpha.lookup(s);
      ++col;
      if(col>=60) {col=0; os<<endl;}
    }
    if(col>0) os<<endl;
  }

  // Write continuous tracks
  n=schema.getNumContinuous();
  for(int i=0 ; i<n ; ++i) {
    const String &name=schema.getContinuousName(i);
    os<<"%"<<name<<endl;
    Vector<Emission>::const_iterator cur=S.begin(), end=S.end();
    for(; cur!=end ; ++cur) {
      const Emission &e=*cur;
      if(e.getContinuous().getDim()<n) INTERNAL_ERROR;
      os<<e.getContinuous()[i]<<endl;
    }    
  }
}



void EmissionSequence::append(const EmissionSequence &other)
{
  S.append(other.S);
}



void EmissionSequence::scaleAllTracks(double multiplier)
{
  const int n=S.size();
  for(int i=0 ; i<n ; ++i) {
    GSL::Vector &c=S[i].getContinuous();
    c.scale(multiplier);
  }
}



void EmissionSequence::recode(int N)
{
  int numDiscrete=schema.getNumDiscrete();
  for(int i=0 ; i<numDiscrete ; ++i) recode(i,N);
}



void EmissionSequence::recode(int trackNum,int N)
{
  //  if(schema.getAlphabet(trackNum).size()<1) {cout<<"tracknum="<<trackNum<<endl;INTERNAL_ERROR;}

  HigherOrderAlphabet H(schema.getAlphabet(trackNum),N);
  int L=length();
  Sequence seq(Symbol(0),L);
  for(int i=0 ; i<L ; ++i) seq[i]=S[i].getDiscrete(trackNum);
  int minNL=min(N,L);
  for(int i=0 ; i<minNL ; ++i) 
    S[i].getDiscrete(trackNum)=H.lookup(seq,0,i+1);
  for(int i=N ; i<L ; ++i) {
    S[i].getDiscrete(trackNum)=H.lookup(seq,i-N+1,N);
  }
}



void EmissionSequence::unencode(int N)
{
  int L=length();
  int numDiscrete=schema.getNumDiscrete();
  for(int i=0 ; i<numDiscrete ; ++i) {
    HigherOrderAlphabet H(schema.getAlphabet(i),N);
    for(int pos=0 ; pos<L ; ++pos) {
      NmerSymbol &s=S[pos].getDiscrete(i);
      s=H.lastBase(s);
    }	
  }  
}



void EmissionSequence::setFilename(const String &f)
{
  filename=f;
}



const String &EmissionSequence::getFilename() const
{
  return filename;
}


void EmissionSequence::copySeqAndSchema(const EmissionSequence &other)
{
  S=other.S;
  schema=other.schema;
}



EmissionSequence *EmissionSequence::getSubsequence(int begin,int length)
{
  EmissionSequence *other=new EmissionSequence(schema);
  other->filename=filename;
  for(int i=begin ; i<begin+length ; ++i)
    other->S.push_back(S[i]);
  return other;
}



Sequence *EmissionSequence::getDiscreteSeq(int trackID)
{
  const int L=length();
  Sequence *seq=new Sequence(Symbol(0),L);
  for(int i=0 ; i<L ; ++i)
    (*seq)[i]=S[i].getDiscrete(0);
  return seq;
}




void EmissionSequence::getExtrema(int trackID,float &minVal,float &maxVal)
{
  const int L=length();
  minVal=POSITIVE_INFINITY;
  maxVal=NEGATIVE_INFINITY;
  for(int i=0 ; i<L ; ++i) {
    double y=S[i][trackID];
    if(y<minVal) minVal=y;
    if(y>maxVal) maxVal=y;
  }
}



