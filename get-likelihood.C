/**************************************************************
 get-likelihood.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#include <string>
#include <iostream>
#include <fstream>
#include "BOOM/Random.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Vector.H"
#include "BOOM/SparseGraph.H"
#include "BOOM/Constants.H"
#include "BOOM/SumLogProbs.H"
#include "BOOM/Regex.H"
#include "HMM.H"
#include "HMMGraph.H"
#include "ForwardAlgorithm.H"
#include "BackwardAlgorithm.H"
#include "EmissionLoader.H"
using namespace std;

class Application {
  Regex fastbRegex;
public:
  Application() : fastbRegex(".*\\.fastb") {}
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



int Application::go(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"p:");
  if(cmd.numArgs()!=2) 
    throw "\n\
get-likelihood [options] <input.hmm> <diretory-or-fastb-file>\n\
  where:  -p \"0 1 3 2 5 6 4\" = apply specified state mapping\n\
\n\
";
  String infile=cmd.arg(0);
  String dir=cmd.arg(1);
  bool specifiedMapping=cmd.option('p');
  cout.precision(30);
  cerr.precision(30);

  // Misc. initialization
  BOOM::catchFloatOverflow();

  // Load HMM
  HMM hmm=(infile);
  HMMGraph hmmGraph(hmm);
  const int numStates=hmm.countStates();

  // Perform mapping, if desired
  Array1D<int> mapping(numStates);
  if(specifiedMapping) {
    String mapString=cmd.optParm('p');
    mapString.trimWhitespace();
    BOOM::Vector<String> *fields=mapString.getFields();
    if(fields->size()!=numStates) 
      throw "Mapping string contains wrong number of states";
    for(int i=0 ; i<numStates ; ++i) mapping[i]=(*fields)[i].asInt();
    delete fields;

    // Apply the mapping
    HMM newHMM(numStates);
    newHMM.getSchema()=hmm.getSchema();
    for(int q=1 ; q<numStates ; ++q)
      newHMM.getEmissionDistr(q)=hmm.getEmissionDistr(mapping[q]);
    for(int q=0 ; q<numStates ; ++q)
      for(int r=0 ; r<numStates ; ++r)
	newHMM.setTransitionProb(q,r,
          hmm.getTransitionProb(mapping[q],mapping[r]));
    hmm=newHMM;
  }

  // Process the training set
  EmissionLoader loader(hmm.getSchema());
  double LL=0.0;
  if(fastbRegex.match(dir)) { // It's actually a filename, not a directory
    EmissionSequence *S=loader.load(dir);
    S->recode(hmm.getOrder()+1); // ###
    //ForwardAlgorithm F(hmmGraph,*S);
    BackwardAlgorithm B(hmmGraph,*S); // ###
    //LL+=F.getLogP();
    LL+=B.getLogP();
    //cout<<F.getLogP()<<" vs. "<<B.getLogP()<<endl; // ###
  }
  else {
    BOOM::Vector<String> files;
    File::getFileList(dir,files);
    int numFiles=files.size();
    for(int i=0 ; i<numFiles ; ++i) {
      String filename=files[i];
      int L=filename.length();
      if(L<7) continue;
      String extension=filename.substr(L-5);
      extension.toupper();
      if(extension!="FASTB") continue;
      String base=filename.substr(0,L-6);
      //String pathFile=base+".path";
      EmissionSequence *S=loader.load(dir+"/"+filename);
      S->recode(hmm.getOrder()+1); // ###
      //BOOM::Vector<int> &path=*loader.loadPathFile(dir+"/"+pathFile);
      ForwardAlgorithm F(hmmGraph,*S);
      //cout<<F<<endl;
      LL+=F.getLogP();
      delete S;
      //delete &path;
    }
  }
  cout<<LL<<endl;
  
  return 0;
}



