/****************************************************************
 deterministic-trainer.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/

#include "BOOM/String.H"
#include <iostream>
#include <fstream>
#include "BOOM/Random.H"
#include "BOOM/CommandLine.H"
#include "BOOM/DnaAlphabet.H"
#include "BOOM/Sequence.H"
#include "HMM.H"
#include "HMMbuilder.H"
#include "HMMreader.H"
#include "BOOM/FastaReader.H"
#include "BOOM/Map.H"
using namespace std;

class Application
{
  DnaAlphabet alphabet;
  HMM *hmm;
  int numStates, numSymbols;
  Vector<Sequence*> *trainingSet;
  Map<int,Map<int,int> > trans;   // state x state x count
  Map<int,Map<Symbol,int> > emit; // state x symbol x count

  Vector<Sequence*> *load(String);
  void train();
  int findSuccessorState(int currentState,const Symbol &);
public:
  int go(int argc,char *argv[]);
};


int main(int argc,char *argv[])
{
  try
    {
      Application app;
      return app.go(argc,argv);
    }
  catch(const char *msg)
    {
      cerr << msg << endl;
      return -1;
    }
  catch(const String &s)
    {
      cerr << s.c_str() << endl;
      return -1;
    }
  catch(string s)
    {
      cerr << s.c_str() << endl;
      return -1;
    }
  catch(...)
    {
      cerr << "unknown exception caught in main()" << endl;
      return -1;
    }
  return 0;
}



void usage()
{
  throw "deterministic-trainer <structure.hmms> <examples.fasta> <outfile>";
}


int Application::go(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=3) usage();
  String structureFile=cmd.arg(0);
  String trainFilename=cmd.arg(1);
  String outfile=cmd.arg(2);

  randomize();
  numSymbols=alphabet.getNumElements();

  // Load HMM structure
  cerr << "Loading HMM structure..." << endl;
  HMMreader reader(alphabet);
  hmm=reader.read(structureFile);
  numStates=hmm->countStates();

  // Load the training set
  trainingSet=load(trainFilename);

  // Train the HMM
  cerr << "Training..." << endl;
  train();
  hmm->normalizeTransitions();
  hmm->normalizeEmissions();

  // Save the results
  cerr << "Saving HMM..." << endl;
  hmm->save(outfile);
  return 0;
}


Vector<Sequence*> *Application::load(String filename)
{
  Vector<Sequence*> *v=new Vector<Sequence*>;

  FastaReader reader(filename);
  String defline, sequence;
  while(reader.nextSequence(defline,sequence))
    {
      Sequence *seq=new Sequence(sequence,alphabet);
      v->push_back(seq);
    }
  return v;
}



void Application::train()
{
  // Initialize pseudocounts
  for(int i=0 ; i<numStates ; ++i)
    for(int j=0 ; j<numStates ; ++j)
      trans[i][j]=hmm->doesTransitionExist(i,j) ? 1 : 0;

  // Iterate through the training strings to obtain
  // evidence for transitions and emissions
  Vector<Sequence*>::iterator cur=trainingSet->begin(),
    end=trainingSet->end();
  for(; cur!=end ; ++cur)
    {
      Sequence &seq=**cur;
      int seqLen=seq.getLength();
      if(seqLen==0) throw "zero-length string!";
      int currentState=0;
      for(int i=0 ; i<seqLen ; ++i)
	{
	  const Symbol &s=seq[i];
	  int nextState=findSuccessorState(currentState,s);
	  if(!trans[currentState].isDefined(nextState))//trans[currentState].find(nextState)==trans.end())
	    trans[currentState][nextState]=0;
	  ++trans[currentState][nextState];
	  if(!emit[nextState].isDefined(s))//emit[nextState].find(s)==emit.end())
	    emit[nextState][s]=0;
	  ++emit[nextState][s];
	  currentState=nextState;
	}
      if(!trans[currentState].isDefined(0))//trans[currentState].find(0)==trans.end())
	trans[currentState][0]=0;
      ++trans[currentState][0];
      if(!hmm->doesTransitionExist(currentState,0))
	throw String("Evidence for transition from ")+currentState
	  +" to zero, not allowed in structure file";
      if(currentState==0) throw "bad!";
    }

  // Set the transition & emission probabilities based on
  // the counts from the training data
  for(int from=0 ; from<numStates ; ++from)
    {
      if(from>0)
	{
	  int totalEmit=0;
	  for(int i=0 ; i<numSymbols ; ++i)
	    {
	      Symbol s=i;
	      if(emit[from].isDefined(s))//emit[from].find(s)!=emit.end())
		totalEmit+=emit[from][s];
	    }
	  if(totalEmit==0) 
	    {
	      cerr << String("Warning: state ")+from+
		" did not emit any symbols in the training data;"
		"\nthis state may be unreachable" << endl;
	      totalEmit=1;
	    }
	  for(int i=0 ; i<numSymbols ; ++i)
	    {
	      Symbol s=i;
	      hmm->setEmissionProb(from,s,emit[from][s]/float(totalEmit));
	    }
	}
      int totalTrans=0;
      for(int to=0 ; to<numStates ; ++to)
	{
	  if(trans[from].isDefined(to))//trans[from].find(to)!=trans.end())
	    totalTrans+=trans[from][to];
	}
      for(int to=0 ; to<numStates ; ++to)
	if(trans[from].isDefined(to))//trans[from].find(to)!=trans.end())
	  hmm->setTransitionProb(from,to,trans[from][to]/float(totalTrans));
    }
}



int Application::findSuccessorState(int currentState,const Symbol &s)
{
  for(int i=1 ; i<numStates ; ++i)
    if(hmm->doesTransitionExist(currentState,i) && 
       hmm->getEmissionProb(i,s)>0) return i;

  for(int i=1 ; i<numStates ; ++i)
    cout<<"trans "<<currentState<<"->"<<i<<"="
	<<hmm->doesTransitionExist(currentState,i)
	<<", emit("<<alphabet.lookup(s)<<")="
	<<(hmm->getEmissionProb(i,s)>0)<<endl;


  throw String("No successor state of ")+currentState+
		   " can emit symbol \""+alphabet.lookup(s)+"\"";
}

