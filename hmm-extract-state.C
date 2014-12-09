/**************************************************************
 hmm-extract-state.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#include <string>
#include <iostream>
#include <fstream>
#include "BOOM/CommandLine.H"
#include "BOOM/Constants.H"
#include "BOOM/HigherOrderAlphabet.H"
#include "HMM.H"
using namespace std;
using namespace BOOM;

class Application {
  void extractState(HMM &oldHMM,int state,const String &outfile);
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


const char *usage="\n\
hmm-extract-state <in.hmm> <state#> <out.hmm>\n\
  use \"all\" for state# to extract all nonzero states\n\
    (give filestem as third argument)\n\
\n\
";

int Application::go(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  int numArgs=cmd.numArgs();
  if(numArgs!=3) throw usage;
  String inFile=cmd.arg(0);
  String stateStr=cmd.arg(1);
  String outFile=cmd.arg(2);
  cout.precision(4);
  HMM oldHMM(inFile);

  stateStr.toupper();
  if(stateStr=="ALL") {
    int numStates=oldHMM.countStates();
    for(int q=1 ; q<numStates ; ++q)
      extractState(oldHMM,q,outFile+q+".hmm");
  }
  else extractState(oldHMM,stateStr.asInt(),outFile);
}



void Application::extractState(HMM &oldHMM,int stateNum,const String &outFile)
{
  HMM newHMM(2);
  Schema schema=newHMM.getSchema()=oldHMM.getSchema();
  newHMM.setTransitionProb(0,1,1);
  newHMM.setTransitionProb(1,0,1);
  newHMM.setTransitionProb(1,1,0);
  newHMM.setTransitionProb(0,0,0);
  Array2D< Array1D<double> > &newEmit=newHMM.getDiscreteEmitDistr();
  Array2D< Array1D<double> > &oldEmit=oldHMM.getDiscreteEmitDistr();
  int numDiscrete=schema.getNumDiscrete();
  int order=oldHMM.getOrder();
  newHMM.changeOrder(order);
  newEmit.resize(2,numDiscrete);
  for(int i=0 ; i<numDiscrete ; ++i) {
    Alphabet &alpha=schema.getAlphabet(i);
    HigherOrderAlphabet H(alpha,order+1);
    int numNmers=H.getNumNmers();
    newEmit[1][i].resize(numNmers);
    for(NmerSymbol nmer=0 ; nmer<numNmers ; ++nmer) 
      newEmit[1][i][nmer]=oldEmit[stateNum][i][nmer];
  }
  newHMM.copyEmissionDistr(oldHMM,stateNum,1);
  newHMM.save(outFile);
  return 0;
}




