/**************************************************************
 HMMreader.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#include "UnscaledHMMreader.H"
#include "BOOM/File.H"
#include <iostream>
#include "UnscaledHMM.H"
using namespace std;


HMMreader::HMMreader(Alphabet &alphabet)
  : alphabet(alphabet),
    transRegex("(\\d+)\\s*->\\s*(\\d+)\\s*:\\s*(\\S+)"),
    emitRegex("state\\s*(\\d+)\\s*:\\s*(\\S.*\\S)"),
    assignRegex("(.)=(\\S+)")
{
  // ctor
}



HMM *HMMreader::read(const String &filename)
{
  Vector<Transition*> transitions;
  Vector<Emission*> emissions;

  int maxState=0;
  File file(filename);
  while(!file.eof())
    {
      String line=file.readLine();
      line.trimWhitespace();
      if(transRegex.search(line))
	{
	  int from=transRegex[1].asInt(), to=transRegex[2].asInt();
	  float P=transRegex[3].asFloat();
	  if(from>maxState) maxState=from;
	  if(to>maxState) maxState=to;
	  transitions.push_back(new Transition(from,to,P));
	}
      else if(emitRegex.search(line))
	{
	  int state=emitRegex[1].asInt();
	  String distrString=emitRegex[2];
	  Vector<String> *fields=distrString.getFields();
	  Vector<String>::iterator cur=fields->begin();
	  for( ; cur!=fields->end() ; ++cur)
	    {
	      if(assignRegex.search(*cur))
		{
		  char symbol=assignRegex[1][0];
		  float P=assignRegex[2].asFloat();
		  emissions.push_back(new Emission(state,symbol,P));
                  cout<<state<<"->"<<symbol<<"="<<P<<endl;
		}
	      else 
		{
		  (*cur).trimWhitespace();
		  if((*cur).length()>0)
		    throw String("Syntax error in HMM structure file: ")
		      + *cur;
		}
	    }
	  delete fields;
	}
      else if(line.length()>0)
	throw String("Syntax error in HMM structure file: ")+line;
    }
  
  HMM *hmm=new HMM(alphabet,maxState+1);
  Vector<Transition*>::iterator cur=transitions.begin();
  for(; cur!=transitions.end() ; ++cur)
    {
      Transition *trans=*cur;
      hmm->setTransitionProb(trans->from,trans->to,trans->P);
      delete trans;
    }

  if(emissions.size()>0)
    {
      Vector<Emission*>::iterator ecur=emissions.begin();
      for(; ecur!=emissions.end() ; ++ecur)
	{
	  Emission *emission=*ecur;
	  hmm->setEmissionProb(emission->state,
			       alphabet.lookup(emission->symbol),
			       emission->P);
	  delete emission;
	}
    }
  else
    {
      //cerr << "using uniform emission distribution" << endl;
      int alphabetSize=alphabet.getNumElements();
      double p=1.0/alphabetSize;
      for(int i=0 ; i<=maxState ; ++i)
	for(int j=0 ; j<alphabetSize ; ++j)
	  hmm->setEmissionProb(i,j,p);
    }
 
  return hmm;
}



