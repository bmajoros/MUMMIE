/*
 scan-with-chains.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
*/
#include <math.h>
#include <iostream>
#include "BOOM/CommandLine.H"
#include "BOOM/Vector.H"
#include "BOOM/Sequence.H"
#include "BOOM/Array2D.H"
#include "BOOM/Regex.H"
#include "BOOM/File.H"
#include "BOOM/SummaryStats.H"
#include "FastViterbi.H"
#include "EmissionLoader.H"
using namespace std;


class Application
{
  Regex deflineRegex;
public:
  Application();
  int main(int argc,char *argv[]);
};


int main(int argc,char *argv[])
  {
    try
      {
	Application app;
	return app.main(argc,argv);
      }
    catch(const char *p)
      {
	cerr << p << endl;
      }
    catch(const string &msg)
      {
	cerr << msg.c_str() << endl;
      }
    catch(const exception &e)
      {
	cerr << "STL exception caught in main:\n" << e.what() << endl;
      }
    catch(...)
      {
	cerr << "Unknown exception caught in main" << endl;
      }
    return -1;
  }



Application::Application()
  : deflineRegex("\\s*>(\\S+)")
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"pg:P");
  if(cmd.numArgs()!=6)
    throw string(
"scan-with-chains <*.hmm> <fg-state> <bg-state> <threshold> <window-size> <*.fastb>"
);
  String hmmFile=cmd.arg(0);
  int fgState=cmd.arg(1).asInt();
  int bgState=cmd.arg(2).asInt();
  float threshold=cmd.arg(3).asFloat();
  int windowSize=cmd.arg(4).asInt();
  String seqFile=cmd.arg(5);

  // Load hmm
  HMM hmm(hmmFile);
  Schema schema=hmm.getSchema();
  hmm.dropContinuousTracks();
  int numStates=hmm.countStates();

  // Load sequence
  EmissionLoader loader(schema);
  EmissionSequence &seq=*loader.load(seqFile);
  const int seqLen=seq.length();
  Regex r("([^\/]+)\\.fastb");
  if(!r.search(seqFile)) throw String("must be a fastb file: ")+seqFile;
  String substrate=r[1];
  //seq.recode(hmm.getOrder()+1);
  //Sequence &S=*seq.getDiscreteSeq(0);

  // Scan
  int lastPos=seqLen-windowSize;
  for(int pos=0 ; pos<=lastPos ; ++pos) {
    EmissionSequence *subseq=seq.getSubsequence(pos,windowSize);
    subseq->recode(hmm.getOrder()+1);
    double fgScore=hmm.getEmissionProb(fgState,*subseq,true);
    double bgScore=hmm.getEmissionProb(bgState,*subseq,true);
    double score=fgScore-bgScore;
    if(score>=threshold) 
      cout<<substrate<<"\tbinding\tsite\t"<<pos+1<<"\t"<<pos+windowSize
	  <<"\t"<<score<<"\t+\t0"<<endl;
    delete subseq;
  }

  return 0;
}



