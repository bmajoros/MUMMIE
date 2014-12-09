/****************************************************************
 TieProfile.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "TieProfile.H"
#include "BOOM/File.H"
using namespace std;
using namespace BOOM;

TieProfile::TieProfile()
{
  // ctor
}



TieProfile::TieProfile(const String &filename)
{
  load(filename);
}



int TieProfile::numEvents() const
{
  return events.size();
}



const TieEvent &TieProfile::operator[](int i) const
{
  return events[i];
}



void TieProfile::load(const String &filename)
{
  File f(filename);
  while(!f.eof()) {
    String line=f.getline();
    if(line.isEmpty()) continue;
    line.trimWhitespace();
    if(line.length()>0 && line[0]=='#') continue;
    line.tolower();
    BOOM::Vector<String> &fields=*line.getFields(" \t,->");
    const int numFields=fields.size();
    if(numFields<2) throw "syntax error in tie profile";
    TieOp op=stringToTieOp(fields[0]);
    if(op==TIE_PARMS) parseTie(line,fields);
    else parseFix(line,fields);
    delete &fields;
  }
}



void TieProfile::parseTie(const String &line,
			  const BOOM::Vector<String> &fields)
{
  const int numFields=fields.size();
  TieEntity entity=stringToTieEntity(fields[1]);
  if(entity==TIE_TRANSITIONS) {
    if(numFields%2) error(line);
    TieTransitionEvent event(TIE_PARMS);
    for(int i=2 ; i<numFields ; i+=2) {
      STATE from=fields[i].asInt(), to=fields[i+1].asInt();
      event.addTransition(from,to);
    }
    events.push_back(event);
  }
  else {
    if(fields[2]!="between") 
      error(line,"expecting keyword \"between\"");
    const TieLevel level=stringToTieLevel(fields[3]);
    TieEvent event(TIE_PARMS,entity);
    if(level==TIE_STATES) {
      if(!canTieStates(entity)) 
	error(line,"can't tie that between states");
      for(int i=4 ; i<numFields ; ++i)
	event.addState(fields[i].asInt());
    }
    else if(!canTieComponents(entity))
      error(line,"can't tie that between components");
    events.push_back(event);
  }  
}



void TieProfile::parseFix(const String &line,
			  const BOOM::Vector<String> &fields)
{
  const int numFields=fields.size();
  TieEntity entity=stringToTieEntity(fields[1]);
  if(entity==TIE_TRANSITIONS) {
    if(numFields%2) error(line);
    TieTransitionEvent event(FIX_PARMS);
    for(int i=2 ; i<numFields ; i+=2) {
      STATE from=fields[i].asInt(), to=fields[i+1].asInt();
      event.addTransition(from,to);
    }
    if(numFields==2) event.addTransition(-1,-1); // means "all transitions"
    events.push_back(event);
    return;
  }
  if(numFields==2 && (entity==TIE_MEANS || 
		      entity==TIE_COV_MATRIX ||
		      entity==TIE_MIXTURE_WEIGHTS || 
		      entity==TIE_MARKOV_CHAINS)) {
    events.push_back(TieEvent(FIX_PARMS,entity));
    return;
  }
  if(entity==TIE_MEANS || entity==TIE_COV_MATRIX)
    error(line,"fix statement has too many options");
  if(numFields<5) error(line);
  if(fields[2]!="in") error(line,"expecting keyword \"in\"");
  const TieLevel level=stringToTieLevel(fields[3]);
  if(level!=TIE_STATES) error(line);
  TieEvent event(FIX_PARMS,entity);
  if(!canFixInStates(entity)) error(line,"can't fix that in states");
  for(int i=4 ; i<numFields ; ++i) event.addState(fields[i].asInt());
  events.push_back(event);
}



void TieProfile::error(const String &line,const String &extra)
{
  throw String("Syntax error in tie profile:\n >> ")+line+"\n"+extra;
}






