/*
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
*/
#include <string>
#include <iostream>
#include <fstream>
#include "BOOM/Random.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Vector.H"
#include "BOOM/SparseGraph.H"
#include "BOOM/Constants.H"
#include "BOOM/SummaryStats.H"
#include "BOOM/Regex.H"
#include "BaumWelch.H"
#include "BaumWelchMT.H"
#include "EmissionLoader.H"
#include "SequenceSet.H"
#include "TieProfile.H"
using namespace std;

const int WINDOW_SIZE=100;

class Application
{
  Alphabet alphabet;
  SequenceSet trainingSet;
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



int Application::go(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=5) 
    throw "\n\
score-groups <training-dir> <hmm> <fg-state> <bg-state> <signal-track-name>\n\
";
  String dir=cmd.arg(0);
  String hmmFile=cmd.arg(1);
  int fgState=cmd.arg(2).asInt();
  int bgState=cmd.arg(3).asInt();
  String sigTrackName=cmd.arg(4);

  // Misc. initialization
  cout.precision(10);
  cerr.precision(10);
  BOOM::catchFloatOverflow();
  Regex r("([^\/]+)\\.fastb");

  // Load hmm
  HMM hmm(hmmFile);
  Schema &schema=hmm.getSchema();
  alphabet=schema.getAlphabet(0);
  const int sigID=schema.lookupContinuousID(sigTrackName);

  // Load the training set
  cerr<<"loading training data"<<endl;
  trainingSet.load(dir,schema);
  trainingSet.recode(hmm.getOrder()+1);

  // Process sequences
  cerr<<"processing sequences"<<endl;
  BOOM::Vector<float> all, cluster;
  const int N=trainingSet.size();
  int nextSeqID=1;
  for(int i=0 ; i<N ; ++i) {
    EmissionSequence &S=*trainingSet[i];
    if(!r.search(S.getFilename())) throw S.getFilename();
    String substrate=r[1];
    int begin;
    double prevX=0;
    const int L=S.length();
    double sum=0;
    int totalLen=0;
    double maxX=NEGATIVE_INFINITY;
    int maxXpos=-1;
    for(int pos=0 ; pos<L ; ++pos) {
      const Emission &e=S[pos];
      const double x=e[sigID];
      if(prevX==0 && x>0) { begin=pos; maxX=x; maxXpos=pos; }
      else if(prevX>0 && x>maxX ) { maxX=x; maxXpos=pos; }
      else if(prevX>0 && x==0) {
	//int start=begin, stop=pos;
	int start=maxX-10-5, stop=maxX+10-5;
	if(start<0) start=0;
	if(stop>L) stop=L;
	for(int j=start ; j<stop ; ++j) {
	  double fgScore=hmm.getEmissionProb(fgState,e,true);
	  double bgScore=hmm.getEmissionProb(bgState,e,true);
	  sum+=fgScore-bgScore;
	  ++totalLen;
	}
      }
      prevX=x;
    }
    sum/=totalLen;
    cout<<sum<<"\t"<<substrate<<endl;
  }
  return 0;
}


