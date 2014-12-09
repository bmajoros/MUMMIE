/**************************************************************
 fastb-length.C
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
#include "HMM.H"
#include "HMMGraph.H"
#include "ForwardAlgorithm.H"
#include "BackwardAlgorithm.H"
#include "EmissionLoader.H"
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
  CommandLine cmd(argc,argv,"d");
  if(cmd.numArgs()!=2) 
    throw "\n\
fastb-length [options] <*.fastb> <*.schema>\n\
  where:  -d = filename is actually a directory\n\
\n\
";
  String infile=cmd.arg(0);
  String schemaFile=cmd.arg(1);
  bool isDir=cmd.option('d');
  cout.precision(30);
  cerr.precision(30);

  // Misc. initialization
  BOOM::catchFloatOverflow();

  // Process the fastb file(s)
  Schema schema(schemaFile);
  EmissionLoader loader(schema);
  if(isDir) {
    const String dir=infile;
    BOOM::Vector<String> files;
    File::getFileList(dir,files);
    int numFiles=files.size();
    for(int i=0 ; i<numFiles ; ++i) {
      String filename=files[i];
      cout<<filename<<" ";
      int L=filename.length();
      if(L<7) continue;
      String extension=filename.substr(L-5);
      extension.toupper();
      if(extension!="FASTB") continue;
      EmissionSequence *S=loader.load(dir+"/"+filename);
      int len=S->length();
      delete S;
      cout<<len<<endl;
    }
  }
  else {
    EmissionSequence *S=loader.load(infile);
    int len=S->length();
    delete S;
    cout<<len<<endl;
  }
  
  return 0;
}



