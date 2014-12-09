/*
 get-extrema.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
*/
#include <math.h>
#include <iostream>
#include "BOOM/CommandLine.H"
#include "BOOM/Vector.H"
#include "BOOM/Sequence.H"
#include "BOOM/Regex.H"
#include "BOOM/File.H"
#include "BOOM/SummaryStats.H"
#include "BOOM/Constants.H"
#include "BOOM/DnaAlphabet.H"
#include "EmissionLoader.H"
using namespace std;

BOOM::Alphabet alphabet;

class Application {
  Regex deflineRegex;
  int trackID;
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
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=3)
    throw string("get-extrema <fastb-dir> <schema.txt> <track-name>");
  String dir=cmd.arg(0);
  String schemaFile=cmd.arg(1);
  String trackName=cmd.arg(2);

  Schema schema(schemaFile);
  EmissionLoader loader(schema);
  trackID=schema.lookupContinuousID(trackName);

  // Iterate through files
  BOOM::Vector<String> files;
  File::getFileList(dir,files);
  const int numFiles=files.size();
  Regex r("([^\/]+)\\.fastb");
  float minVal=POSITIVE_INFINITY, maxVal=NEGATIVE_INFINITY;
  for(int i=0 ; i<numFiles ; ++i) {
    String seqFile=files[i];
    if(!r.search(seqFile)) continue;

    // Load sequence
    String fullPath=dir+"/"+seqFile;
    EmissionSequence &seq=*loader.load(fullPath);
    const int seqLen=seq.length();

    // Get stats
    float thisMin, thisMax;
    seq.getExtrema(trackID,thisMin,thisMax);
    if(thisMin<minVal) minVal=thisMin;
    if(thisMax>maxVal) maxVal=thisMax;

    delete &seq;
  }
  cout<<minVal<<" - "<<maxVal<<endl;

  return 0;
}



