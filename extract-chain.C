/**************************************************************
 extract-chain.C : For training an HMM using Expectation Maximization
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#include <string>
#include <iostream>
#include <fstream>
#include "BOOM/CommandLine.H"
#include "BOOM/VectorSorter.H"
#include "BOOM/Constants.H"
#include "BOOM/HigherOrderAlphabet.H"
#include "BOOM/NthOrderStringIterator.H"
#include "EmissionLoader.H"
#include "SequenceSet.H"
#include "HMM.H"
using namespace std;


class Application {
public:
  int go(int argc,char *argv[]);
protected:
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
  CommandLine cmd(argc,argv,"as:L:i");
  if(cmd.numArgs()!=2)
    throw "\n\
extract-chain <*.hmm> <state>\n\
\n\
";
  String hmmFile=cmd.arg(0);
  int state=cmd.arg(1).asInt();

  // Misc. initialization
  cout.precision(10);
  cerr.precision(10);
  BOOM::catchFloatOverflow();

  // Load HMM
  HMM hmm(hmmFile);
  Schema &schema=hmm.getSchema();
  const int order=hmm.getOrder();

  // Extract chain
  Array2D< Array1D<double> > &E=hmm.getDiscreteEmitDistr();
  int numDiscrete=schema.getNumDiscrete();
  for(int i=0 ; i<numDiscrete ; ++i) {
    cout<<order<<" order"<<endl;
    Alphabet &alphabet=schema.getAlphabet(i);
    cout<<"alphabet:\n";
    alphabet.save(cout);
    HigherOrderAlphabet H(alphabet,order+1);
    const int numNmers=H.getNumNmers();
    for(NmerSymbol nmer=1 ; nmer<numNmers ; ++nmer) {
      Sequence S=H.lookup(nmer);
      S.printOn(cout,alphabet);
      cout<<"\t"<<E[state][i][nmer]<<endl;
    }
  }	

  return 0;
}



