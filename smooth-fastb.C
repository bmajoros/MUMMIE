/****************************************************************
 smooth-fastb.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/FastbReader.H"
#include "BOOM/Exceptions.H"
#include "BOOM/Set.H"
using namespace std;
using namespace BOOM;


class Application {
  Set<String> tracks;
  bool notAllTracks;
  void smooth(FastbContinuousSeq &,int windowSize,int numIterations);
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
  CommandLine cmd(argc,argv,"t:");
  if(cmd.numArgs()!=4)
    throw String("\n\
smooth-fastb [ops] <in.fastb> <window-size> <#iterations> <out.fastb>\n\
    where: -t abc,xyz,... = smooth only the listed tracks\n\
");
  const String infile=cmd.arg(0);
  const int windowSize=cmd.arg(1).asInt();
  const int numIterations=cmd.arg(2).asInt();
  const String outfile=cmd.arg(3);
  notAllTracks=cmd.option('t');
  if(notAllTracks) {
    String trackList=cmd.optParm('t');
    BOOM::Vector<String> &trackV=*trackList.getFields(",");
    int n=trackV.size();
    for(int i=0 ; i<n ; ++i) tracks.insert(trackV[i]);
    delete &trackV;
  }

  ofstream os(outfile.c_str());
  FastbReader reader(infile);
  while(!reader.eof()) {
    FastbSequence *seq=reader.nextSequence();
    if(!seq) break;
    if(seq->isContinuous()) 
      smooth(static_cast<FastbContinuousSeq&>(*seq),windowSize,numIterations);
    os<<*seq;
    delete seq;
  }	

  return 0;
}



void Application::smooth(FastbContinuousSeq &seq,int windowSize,
			 int numIterations)
{
  if(notAllTracks && !tracks.isMember(seq.getID())) return;
  DoubleVector &S=seq.getSeq();
  int L=S.size();
  DoubleVector newSeq=S;
  int halfWindow=windowSize/2;
  for(int i=0 ; i<numIterations ; ++i) {
    for(int pos=1 ; pos<L-1 ; ++pos) {
      double delta=halfWindow;
      if(delta>pos) delta=pos;
      else if(delta>L-1-pos) delta=L-1-pos;
      int first=pos-delta, last=pos+delta;
      double ave=0, d=last-first+1;
      for(int j=first ; j<=last ; ++j) 
	//ave+=(S[pos]+S[j])/(2*d);
	ave+=S[j]/d;
      newSeq[pos]=ave;
    }
    S=newSeq;
  }
}



