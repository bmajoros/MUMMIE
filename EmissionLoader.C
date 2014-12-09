/****************************************************************
 EmissionLoader.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "EmissionLoader.H"
#include "BOOM/FastaReader.H"
#include "BOOM/FastcReader.H"
#include "BOOM/FastbReader.H"
using namespace std;
using namespace BOOM;

EmissionLoader::EmissionLoader(Schema &schema)
  : length(0), schema(schema)
{
  // ctor
}



void EmissionLoader::loadContinuous(const String &filename)
{
  FastcReader cf(filename);
  String defline, id;
  int N=continuous.size();
  continuous.push_back(Vector<double>());
  Vector<double> &V=continuous[N];
  if(!cf.nextSequence(defline,id,V))
    throw String("Error loading ")+filename+" in EmissionLoader";
  continuousNames.push_back(id);
  int s=V.size();
  if(s>length) length=s;
}



void EmissionLoader::loadDiscrete(const String &filename,Alphabet &alpha)
{
  FastaReader reader(filename,alpha);
  String defline, seq, id, remainder;
  if(!reader.nextSequence(defline,seq))
    throw String("Error loading ")+filename+" in EmissionLoader";
  FastaReader::parseDefline(defline,id,remainder);
  int N=discrete.size();
  discrete.push_back(Sequence());
  Sequence &S=discrete[N];
  S.copyFrom(seq,schema.getAlphabet(schema.lookupDiscreteID(id)));
  discreteNames.push_back(id);
  int s=S.getLength();
  if(s>length) length=s;
}



void EmissionLoader::loadMixed(const String &filename)
{
  FastbReader reader(filename);
  FastbSequence *seq;
  while(seq=reader.nextSequence()) {
    if(seq->getType()==FASTB_DISCRETE) {
      String id=seq->getID();
      discreteNames.push_back(id);
      Alphabet &alpha=schema.getAlphabet(schema.lookupDiscreteID(id));
      Sequence S(static_cast<FastbDiscreteSeq*>(seq)->getSeq(),alpha);
      discrete.push_back(S);
      if(S.getLength()>length) length=S.getLength();
    }
    else {
      continuousNames.push_back(seq->getID());
      const DoubleVector &V=static_cast<FastbContinuousSeq*>(seq)->getSeq();
      continuous.push_back(V);
      if(V.size()>length) length=V.size();
    }
    delete seq;
  }
}



EmissionSequence *EmissionLoader::compile()
{
  // Allocate emission object
  int numD=discrete.size(), numC=continuous.size();
  EmissionSequence *E=new EmissionSequence(schema);

  // Check lengths
  for(int i=0 ; i<numD ; ++i) 
    if(discrete[i].getLength()!=length)
      throw discreteNames[i]+" has nonmatching length: "+
	discrete[i].getLength()+" vs. "+length;
  for(int i=0 ; i<numC ; ++i) 
    if(continuous[i].size()!=length)
      throw continuousNames[i]+" has nonmatching length: "+
	continuous[i].size()+" vs. "+length;

  // Create emission points
  int dPerm[numD], cPerm[numC];
  for(int i=0 ; i<numD ; ++i) {
    dPerm[i]=schema.lookupDiscreteID(discreteNames[i]);
    if(dPerm[i]<0) throw String("Track not in schema: ")+discreteNames[i];
  }	
  for(int i=0 ; i<numC ; ++i) {
    cPerm[i]=schema.lookupContinuousID(continuousNames[i]);
    if(cPerm[i]<0) throw String("track not in schema : ")+continuousNames[i];
  }
  for(int pos=0 ; pos<length ; ++pos) {
    Emission e;
    e.resize(numD,numC);
    for(int i=0 ; i<numD ; ++i)
      e.getDiscrete(dPerm[i])=discrete[i][pos];
    GSL::Vector &c=e.getContinuous();
    for(int i=0 ; i<numC ; ++i)
      c[cPerm[i]]=continuous[i][pos];
    E->append(e);
  }

  // Re-init for next time
  discrete.clear();
  continuous.clear();
  discreteNames.clear();
  continuousNames.clear();
  length=0;
  
  return E;
}



EmissionSequence *EmissionLoader::load(const String &fastbFile)
{
  loadMixed(fastbFile);
  return compile();
}



Vector<int> *EmissionLoader::loadPathFile(const String &filename)
{
  ifstream is(filename.c_str());
  if(!is.good()) throw String("Can't open file: ")+filename;
  Vector<int> *v=new Vector<int>;
  while(!is.eof()) {
    int state=-1;
    is>>state;
    if(is.eof() || state<0) break;
    v->push_back(state);
  }
  return v;
}



