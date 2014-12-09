/**************************************************************
 find-peaks.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#include <string>
#include <iostream>
#include <fstream>
#include "BOOM/CommandLine.H"
#include "BOOM/Vector.H"
#include "BOOM/Array2D.H"
#include "BOOM/Constants.H"
#include "BOOM/Regex.H"
#include "BOOM/SummaryStats.H"
#include "EmissionLoader.H"
using namespace std;

class Application {
  Regex substrateRegex;
  int graphWidth, numContinuous, sampleSize;
  Array2D<float> hist; // [track][pos]
  void report(EmissionSequence &,String substrate,int trackID);
  void tabulate(EmissionSequence &S,int mainTrackID);
  bool isPeak(EmissionSequence &S,int length,int begin,int end,
	      int trackID,double val);
public:
  float threshold;
  Application();
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



Application::Application() 
  : substrateRegex("([^/]+)\\.fastb")
{
}



int Application::go(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"dg:t:");
  if(cmd.numArgs()!=3) 
    throw "\n\
find-peaks [options] <*.fastb> <*.schema> <track-name>\n\
  where:  -d = filename is actually a directory\n\
          -t T = apply threshold T (default=0.1)\n\
          -g W = emit histograms around peak (+/- width W)\n\
\n\
";
  String infile=cmd.arg(0);
  String schemaFile=cmd.arg(1);
  String trackName=cmd.arg(2);
  bool isDir=cmd.option('d');
  bool wantGraph=cmd.option('g');
  graphWidth=wantGraph ? cmd.optParm('g').asInt() : 0;
  threshold=cmd.option('t') ? cmd.optParm('t').asFloat() : 0.1;

  // Process the fastb file(s)
  Schema schema(schemaFile);
  numContinuous=schema.getNumContinuous();
  if(wantGraph) { 
    hist.resize(numContinuous,2*graphWidth+1); 
    hist.setAllTo(0); 
    sampleSize=0;
  }
  EmissionLoader loader(schema);
  int trackID=schema.lookupContinuousID(trackName);
  if(isDir) {
    const String dir=infile;
    BOOM::Vector<String> files;
    File::getFileList(dir,files);
    int numFiles=files.size();
    for(int i=0 ; i<numFiles ; ++i) {
      String filename=files[i];
      int L=filename.length();
      if(L<7) continue;
      String extension=filename.substr(L-5);
      extension.toupper();
      if(extension!="FASTB") continue;
      EmissionSequence *S=loader.load(dir+"/"+filename);
      if(!substrateRegex.search(filename)) throw filename+" : can't parse";
      String substrate=substrateRegex[1];
      if(wantGraph) tabulate(*S,trackID);
      else report(*S,substrate,trackID);
      delete S;
    }
    if(wantGraph) {
      int W=2*graphWidth+1;
      for(int j=0 ; j<numContinuous ; ++j) {
	//float sum=0;
	//for(i=0 ; i<W ; ++i) sum+=hist[j][i];
	for(i=0 ; i<W ; ++i) {
	  //hist[j][i]/=sum;
	  hist[j][i]/=sampleSize;
	  int x=i-graphWidth;
	  cout<<x<<"\t"<<hist[j][i]<<endl;
	}
	cout<<endl;
      }
    }
  }
  else {
    EmissionSequence *S=loader.load(infile);
    if(!substrateRegex.search(infile)) throw infile+" : can't parse";
    String substrate=substrateRegex[1];
    report(*S,substrate,trackID);
    delete S;
  }
  
  return 0;
}


void Application::tabulate(EmissionSequence &S,int mainTrackID)
{
  int len=S.length();
  int b=0, e=0;
  while(b<len) {
    double x=S[b][mainTrackID];
    while(e+1<len && S[e+1][mainTrackID]==x) ++e;
    if(isPeak(S,len,b,e,mainTrackID,x) && x>=threshold) {
      int peak=(b+e)/2;
      int left=peak-graphWidth; 
      int right=peak+graphWidth;
      if(left>=0 && right<len) {
	for(int i=left ; i<=right ; ++i) {
	  const Emission &e=S[i];
	  for(int j=0 ; j<numContinuous ; ++j) {
	    float val=e[j];
	    hist[j][i-left]+=val;
	  }
	}
	++sampleSize;
      }
    }
    b=e+1;
    e=b;
  }
}



void Application::report(EmissionSequence &S,String substrate,int trackID)
{
  int len=S.length();
  int b=0, e=0;
  while(b<len) {
    double x=S[b][trackID];
    while(e+1<len && S[e+1][trackID]==x) ++e;
    if(isPeak(S,len,b,e,trackID,x) && x>=threshold) {
      int pos=(b+e)/2;
      cout<<substrate<<"\tPAR\tpeak\t"<<pos+1<<"\t"<<pos+1<<"\t"<<x<<"\t+\t.\n";
    }
    b=e+1;
    e=b;
  }
}



bool Application::isPeak(EmissionSequence &S,int length,int begin,int end,
			 int trackID,double val)
{
  if(begin>0 && S[begin-1][trackID]>val) return false;
  if(end+1<length && S[end+1][trackID]>val) return false;
  return true;
}




