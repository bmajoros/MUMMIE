/****************************************************************
 sample.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Random.H"
#include "BOOM/Constants.H"
#include "BOOM/DnaAlphabet.H"
#include "BOOM/SummaryStats.H"
#include "HMM.H"
using namespace std;
using namespace BOOM;


class Application {
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
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
  {
    // Process command line
    CommandLine cmd(argc,argv,"m:M:r:S:DR:");
    if(cmd.numArgs()!=3)
      throw String(
"\n\
sample [options] <in.hmm> <out.fastb> <out.paths>\n\
   where options are:\n\
      -m <min-length>\n\
      -M <max-length>\n\
      -r <randomization-seed>\n\
      -S X = stay in state X\n\
      -D = discrete tracks only\n\
      -R N = repeat sampling N times\n\
");
    String hmmFile=cmd.arg(0);
    String outfile=cmd.arg(1);
    String pathFile=cmd.arg(2);
    int minLength=cmd.option('m') ? cmd.optParm('m').asInt() : 1;
    int maxLength=cmd.option('M') ? cmd.optParm('M').asInt() : LARGEST_INTEGER;
    bool restrictState=cmd.option('S');
    int wantState=restrictState ? cmd.optParm('S').asInt() : 0;
    bool discreteOnly=cmd.option('D');
    int numSamples=cmd.option('R') ? cmd.optParm('R').asInt() : 1;

    // Load HMM
    HMM hmm(hmmFile);
    hmm.normalizeTransitions();
    const int order=hmm.getOrder();
    //cout<<hmm<<endl;

    // Perform sampling
    if(cmd.option('r')) SeedRandomizer(cmd.optParm('r').asUnsigned());
    else randomize();
    ofstream osFastb(outfile.c_str()), osPath(pathFile.c_str());
    for(int sampleNum=0 ; sampleNum<numSamples ; ++sampleNum) {
      BOOM::Vector<STATE> statePath;
      EmissionSequence *seq=NULL;
      const int MAX_SAMPLES=1000;
      BOOM::Vector<int> lengths;
      for(int i=0 ; i<MAX_SAMPLES ; ++i) {
	if(restrictState) {
	  seq=hmm.sampleFromState(wantState,maxLength);
	  lengths.push_back(seq->length());
	  break;
	}
	else {
	  seq=hmm.sample(&statePath,maxLength);
	  if(!seq) continue;
	  int L=seq->length();
	  lengths.push_back(L);
	  if(L>=minLength && L<=maxLength) break;
	  delete seq;
	  statePath.clear();
	  seq=NULL;
	}
      }
      if(!seq) {
	SummaryStats stats(lengths);
	cout<<"mean length = "<<stats.getMean()<<endl;
	throw "Length constraints not satisfied";
      }
      seq->unencode(order+1); // convert from Nmers to individual bases
      if(discreteOnly) seq->getSchema().setNumContinuous(0);
      seq->save(osFastb);
      delete seq;
      BOOM::Vector<STATE>::iterator cur=statePath.begin(), end=statePath.end();
      for(; cur!=end ; ++cur) osPath<<*cur<<"\n";
      osPath<<endl;
    }
    
    // Clean up
    return 0;
  }

