/**************************************************************
 log-likelihood.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/

#include <string>
#include <iostream>
#include <fstream>
#include "BOOM/Vector.H"
#include "BOOM/Regex.H"
#include "BOOM/File.H"
#include "BOOM/CommandLine.H"
#include "HMM.H"
#include "FastViterbi.H"
using namespace std;
using namespace BOOM;

class Application
{
  BOOM::Vector< BOOM::Vector<double>*> trainingSet;
  Regex deflineRegex;
  void load(String filename);
public:
  Application();
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
    throw "log-likelihood <*.hmm> <datafile>";
}



Application::Application()
  : deflineRegex("\\s*>(\\S+)")
{
  // ctor
}



int Application::go(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=2) usage();
  String hmmFile=cmd.arg(0);
  String dataFile=cmd.arg(1);

  // Load HMM
  HMM hmm(hmmFile);
  FastViterbi viterbi(hmm);
  
  // Process the FASTA file
  load(dataFile);
  int numTrain=trainingSet.size();
  int n=trainingSet.size();
  for(int i=0 ; i<n ; ++i) {
    BOOM::Vector<double> &seq=*trainingSet[i];
    double pathScore;
    BOOM::Vector<int> &path=*viterbi.getPath(seq,pathScore);
    cout<<pathScore<<endl;
  }
  
  return 0;
}


void Application::load(String filename)
{
  File file(filename);
  while(!file.eof()) {
    String line=file.getline();
    if(deflineRegex.match(line)) {
      line=file.getline();
      BOOM::Vector<BOOM::String> *fields=line.getFields();
      int numFields=fields->size();
      BOOM::Vector<double> *seq=new BOOM::Vector<double>;
      for(int i=0 ; i<numFields ; ++i)
	seq->push_back((*fields)[i].asDouble());
      delete fields;
      trainingSet.push_back(seq);
    }
  }
}
