/**************************************************************
 extract-state-sequences.C
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
#include "BOOM/Constants.H"
#include "EmissionLoader.H"
#include "SequenceSet.H"
#include "HMM.H"
using namespace std;

class Application
{
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
extract-state-sequences <in.hmm> <in.fastb> <track-name> <in.path> <state>\n\
";
  String hmmFile=cmd.arg(0);
  String fastbFile=cmd.arg(1);
  String trackName=cmd.arg(2);
  String pathFile=cmd.arg(3);
  const int wantState=cmd.arg(4).asInt();

  // Load HMM
  HMM hmm(hmmFile);
  Schema &schema=hmm.getSchema();

  // Load sequence
  EmissionLoader loader(schema);
  EmissionSequence &S=*loader.load(fastbFile);
  BOOM::Vector<int> &path=*loader.loadPathFile(pathFile);
  const int trackID=schema.lookupDiscreteID(trackName);
  Alphabet alpha=schema.getAlphabet(trackID);
  
  // Process path
  int seqID=1;
  const int L=path.size();
  for(int i=0 ; i<L ; ++i) {
    if(path[i]==wantState) {
      int begin=i;
      for(; i<L ; ++i) if(path[i]!=wantState) break;
      int end=i-1;
      cout<<">"<<seqID++<<endl;
      for(int j=begin ; j<=end ; ++j) {
	cout<<alpha.lookup(S[j].getDiscrete(trackID));
	if(j>begin && (j-begin)%60==0) cout<<endl;
      }
      cout<<endl;
    }
  }

  delete &S;
  delete &path;
  return 0;
}


