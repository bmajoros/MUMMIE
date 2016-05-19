/****************************************************************
 Schema.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "Schema.H"
#include "BOOM/Exceptions.H"
#include "BOOM/File.H"
using namespace std;


Schema::Schema(int numDiscrete,int numContinuous)
  : numDiscrete(numDiscrete), numContinuous(numContinuous),
    discreteNames(numDiscrete), continuousNames(numContinuous),
    alphabets(numDiscrete)
{
  // ctor
}



Schema::Schema(const String &filename)
{
  ifstream is(filename.c_str());
  if(!is.good() || !File::exists(filename)) throw filename+" : can't open";
  load(is);
}



const String &Schema::getDiscreteName(int i)
{
  return discreteNames[i];
}



const String &Schema::getContinuousName(int i)
{
  return continuousNames[i];
}



void Schema::setDiscreteName(int i,const String &name)
{
  discreteNames[i]=name;
  discreteIDs[name]=i;
}



void Schema::setContinuousName(int i,const String &name)
{
  continuousNames[i]=name;
  continuousIDs[name]=i;
}



int Schema::getNumDiscrete() const
{
  return numDiscrete;
}



int Schema::getNumContinuous() const
{
  return numContinuous;
}



void Schema::setNumDiscrete(int d)
{
  numDiscrete=d;
  discreteNames.resize(d);
  alphabets.resize(d);
}



void Schema::setNumContinuous(int c)
{
  numContinuous=c;
  continuousNames.resize(c);
}



void Schema::load(istream &is)
{
  Vector<String> d, c;
  String line;
  while(!is.eof()) {
    line.getline(is);
    line.trimWhitespace();
    Vector<String> &fields=*line.getFields(" \t:=");
    int n=fields.size();
    if(n>1) {
      if(fields[1]=="discrete") {
	if(n!=3) throw String("Error in schema: ")+line;
	const int id=d.size();
	const String name=fields[0];
	d.push_back(name);
	discreteIDs[name]=id;
	alphabets.push_back(Alphabet(fields[2].c_str()));
      }
      else if(fields[1]=="continuous") {
	const int id=c.size();
	const String name=fields[0];
	c.push_back(name);
	continuousIDs[name]=id;
      }
      else throw String("Error in schema file: ")+line;
    }
    else if(n>0) throw String("Error in schema file: ")+line;
    delete &fields;
  }
  numDiscrete=d.size();
  numContinuous=c.size();
  discreteNames=d;
  continuousNames=c;
}



void Schema::save(ostream &os) const
{
  for(int i=0 ; i<numDiscrete ; ++i) {
    os<<discreteNames[i]<<" : discrete ";alphabets[i].save(os);os<<"\n";
  }
  for(int i=0 ; i<numContinuous ; ++i)
    os<<continuousNames[i]<<" : continuous\n";
}



Alphabet &Schema::getAlphabet(int i)
{
  return alphabets[i];
}



int Schema::lookupDiscreteID(const String &name)
{
  if(!discreteIDs.isDefined(name)) return -1;
  return discreteIDs[name];
}


int Schema::lookupContinuousID(const String &name)
{
  if(!continuousIDs.isDefined(name)) return -1;
  return continuousIDs[name];
}



void Schema::addContinuousTrack(const String &name)
{
  continuousNames.push_back(name);
  continuousIDs[name]=numContinuous;
  ++numContinuous;
}



void Schema::dropContinuousTrack(const String &name)
{
  int id=continuousIDs[name];
  continuousNames.cut(id);
  continuousIDs.clear();
  --numContinuous;
  for(int i=0 ; i<numContinuous ; ++i)
    continuousIDs[continuousNames[i]]=i;
}



void Schema::dropDiscreteTrack(const String &name)
{
  int id=continuousIDs[name];
  continuousNames.cut(id);
  alphabets.cut(id);
  continuousIDs.clear();
  --numDiscrete;
  for(int i=0 ; i<numDiscrete ; ++i)
    continuousIDs[continuousNames[i]]=i;
}


