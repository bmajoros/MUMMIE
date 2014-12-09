/**************************************************************
 install-chain.C
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
  Application();
  int go(int argc,char *argv[]);
protected:
  Regex rxOrder, rxAlphabet;
  String find(Regex &,istream &);
};



Application::Application()
  : rxOrder("\\s*(\\d+)\\s+order\\s*"),
    rxAlphabet("\\s*alphabet\\s*:\\s*")
{
}




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
  if(cmd.numArgs()!=4)
    throw "\n\
install-chain <*.chain> <in.hmm> <out.hmm> <state>\n\
\n\
";
  String chainFile=cmd.arg(0);
  String hmmFile=cmd.arg(1);
  String outfile=cmd.arg(2);
  int state=cmd.arg(3).asInt();

  // Misc. initialization
  cout.precision(10);
  cerr.precision(10);
  BOOM::catchFloatOverflow();

  // Load HMM
  HMM hmm(hmmFile);
  Schema &schema=hmm.getSchema();
  const int order=hmm.getOrder();
  int numDiscrete=schema.getNumDiscrete();

  // Load chain
  Array2D< Array1D<double> > &E=hmm.getDiscreteEmitDistr();
  ifstream is(chainFile.c_str());
  if(!is.good()) throw String("Error opening file: ")+chainFile;
  String line;
  for(int i=0 ; i<numDiscrete ; ++i) {
    line=find(rxOrder,is);
    const int order=rxOrder[1];
    find(rxAlphabet,is);
    line.getline(is);
    Alphabet alpha(line.c_str());
    schema.getAlphabet(i)=alpha;
    HigherOrderAlphabet H(alpha,order+1);
    const int numNmers=H.getNumNmers();
    Array1D<double> &row=E[state][i];
    row.resize(numNmers);
    for(NmerSymbol nmer=1 ; nmer<numNmers ; ++nmer) {
      line.getline(is);
      line.trimWhitespace();
      if(is.eof() || line.isEmpty()) break;
      BOOM::Vector<BOOM::String> &fields=*line.getFields();
      Sequence seq(fields[0],alpha);
      NmerSymbol nmer=H.lookup(seq);
      row[nmer]=fields[1].asDouble(); // in log space
      delete &fields;
    }
  }

  // Save
  hmm.save(outfile);

  return 0;
}



String Application::find(Regex &regex,istream &is)
{
  String line;
  while(!is.eof()) {
    line.getline(is);
    if(regex.match(line)) return line;
  }
  throw String("Syntax error in chain file; expecting: ")+regex;
}




