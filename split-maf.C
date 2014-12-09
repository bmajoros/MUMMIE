/*
 split-maf.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
*/
#include <math.h>
#include <fstream>
#include <iostream>
#include "BOOM/CommandLine.H"
#include "BOOM/Exceptions.H"
#include "BOOM/Vector.H"
#include "BOOM/Sequence.H"
#include "BOOM/Array2D.H"
#include "BOOM/Regex.H"
#include "BOOM/File.H"
#include "BOOM/SummaryStats.H"
#include "BOOM/MultiAlignment.H"
using namespace std;

String keep[]={"PanTro2","GorGor1","PonAbe2","RheMac2","PapHam1","CalJac1","TarSyr1","MicMur1","OtoGar1",""};

class Application {
  Regex ENSTregex;
  //  void flush(Vector<MultiAlignment*> &);
  String getTranscriptID(MultiAlignment *maf,int &chunkID);
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
  : ENSTregex("([^_]+)_(\\S+)")
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=2)
    throw string(
"split-maf <in.maf> <outdir>\n\
");
  String inFile=cmd.arg(0);
  String outdir=cmd.arg(1);
  
  int numKeep=0;
  for(int i=0 ; ; ++i, ++numKeep) if(keep[i].isEmpty()) break;

  ifstream is(inFile.c_str());
  MultiAlignment *maf;
  String prevTranscript;
  //Vector<MultiAlignment*> chunks;
  while(maf=MultiAlignment::nextAlignmentFromMAF(is)) {
    maf->toupper();
    int chunkID;
    String transcriptID=getTranscriptID(maf,chunkID);
    //String outfile=outdir+"/"+transcriptID+"_"+chunkID+".maf";
    String outfile=outdir+"/"+transcriptID+".maf";
    ofstream os(outfile.c_str());
    os<<*maf<<endl;
    delete maf;
    /*
    if(transcriptID==prevTranscript) chunks.push_back(maf);
    else {
      flush(chunks);
      chunks.push_back(maf);
    }
    prevTranscript=transcriptID;
    */
  }
  //flush(chunks);

  return 0;
}



String Application::getTranscriptID(MultiAlignment *maf,int &chunkID) 
{
  int numTracks=maf->getNumTracks();
  for(int i=0 ; i<numTracks ; ++i) {
    AlignmentTrack &track=maf->getIthTrack(i);
    const String &name=track.getName();
    if(ENSTregex.match(name)) {
      String transcriptID=ENSTregex[1];
      chunkID=ENSTregex[2].asInt();
      //return transcriptID;
      return name;
    }
  }
  INTERNAL_ERROR;
}



/*
void Application::flush(Vector<MultiAlignment*> &chunks)
{
  int n=chunks.size();
  //cout<<"n="<<n<<endl;
  for(int i=0 ; i<n ; ++i) {
    MultiAlignment *maf=chunks[i];
    int chunkID;
    String transcriptID=getTranscriptID(maf,chunkID);
    cout<<chunkID<<" "<<i<<endl;
    if(chunkID!=i+1) {cout<<chunkID<<" vs "<<i<<endl;throw "bad";}

    delete maf;
  }
  chunks.clear();
}
*/


