/**************************************************************
 HMMreader.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#include "HMMreader.H"
#include "BOOM/File.H"
#include <iostream>
#include "HMM.H"
using namespace std;


HMMreader::HMMreader()
  : transRegex("(\\d+)\\s*->\\s*(\\d+)\\s*:\\s*(\\S+)")
{
  // ctor
}



HMM *HMMreader::read(const String &filename)
{
  BOOM::Vector<Transition*> transitions;
  int maxState=0;
  File file(filename);
  while(!file.eof()) {
    String line=file.readLine();
    line.trimWhitespace();
    if(transRegex.search(line))	{
      int from=transRegex[1].asInt(), to=transRegex[2].asInt();
      float P=transRegex[3].asFloat();
      if(from>maxState) maxState=from;
      if(to>maxState) maxState=to;
      transitions.push_back(new Transition(from,to,P));
    }
    else if(line.length()>0)
      throw String("Syntax error in metamodel file: ")+line;
  }
  
  HMM *hmm=new HMM(maxState+1);
  BOOM::Vector<Transition*>::iterator cur=transitions.begin();
  for(; cur!=transitions.end() ; ++cur) {
    Transition *trans=*cur;
    hmm->setTransitionProb(trans->from,trans->to,trans->P);
    delete trans;
  }

  return hmm;
}



