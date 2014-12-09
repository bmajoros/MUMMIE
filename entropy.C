/**************************************************************
 entropy.C 
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#include <string>
#include <iostream>
#include <fstream>
#include "BOOM/FastaReader.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Constants.H"
#include "BOOM/Exceptions.H"
#include "BOOM/SequenceEntropy.H"
using namespace std;

class Application {
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
  if(cmd.numArgs()!=2)
    throw "\n\
entropy <in.fasta> <order>\n\
\n\
";
  String filename=cmd.arg(0);
  int order=cmd.arg(1).asInt();

  FastaReader reader(filename);
  String def, seq;
  Sequence S;
  while(reader.nextSequence(def,seq)) {
    Sequence s(seq,DnaAlphabet::global());
    S+=s;
  }
  double maxH, maxJH;
  double H=SequenceEntropy::entropy(S,maxH);
  double jH=SequenceEntropy::jointEntropy(S,order,maxJH);
  double cH=SequenceEntropy::conditionalEntropy(S,order);
  cout<<"H="<<H<<" (max="<<maxH<<") jointH="<<jH<<" (max="<<maxJH<<") condH="
      <<cH<<endl;
}
