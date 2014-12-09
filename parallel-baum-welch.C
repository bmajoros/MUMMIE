/**************************************************************
 parallel-baum-welch.C : For training an HMM using Expectation Maximization
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/

#include <string>
#include <iostream>
#include <fstream>
#include "BOOM/FastaReader.H"
#include "BOOM/Sequence.H"
#include "BOOM/Random.H"
#include "BOOM/DnaAlphabet.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Vector.H"
#include "HMMreader.H"
#include "HMMbuilder.H"
#include "ParallelBaumWelch.H"
using namespace std;

class Application
{
  Vector<Sequence*> *trainingSet;
  DnaAlphabet alphabet;

  Vector<Sequence*> *load(String);
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
  catch(...)
    {
      cerr << "unknown exception caught in main()" << endl;
      return -1;
    }
  return 0;
}



void usage()
{
  throw "parallel-baum-welch <structure.hmms> <examples.fasta> <#iterations> <outfile> <num-threads>";
}


int Application::go(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=5) usage();
  String structureFile=cmd.arg(0);
  String trainFilename=cmd.arg(1);
  int maxIterations=cmd.arg(2).asInt();
  String outfile=cmd.arg(3);
  int numThreads=cmd.arg(4).asInt();
  randomize();

  // Load HMM structure
  HMMreader reader(alphabet);
  HMM *hmm=reader.read(structureFile);
  //cerr << "Base model:\n" << *hmm << endl;

  // Load the training set
  trainingSet=load(trainFilename);

  // Train the HMM using the Baum-Welch algorithm
  cerr << "Training..." << endl;
  ParallelBaumWelch(*hmm,maxIterations,*trainingSet,numThreads);
  cerr << "Final model:\n" << *hmm << endl;

  // Save the results
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






