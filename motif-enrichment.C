/*
 motif-enrichment.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
*/
#include <math.h>
#include <iostream>
#include "BOOM/CommandLine.H"
#include "BOOM/Sequence.H"
#include "BOOM/Constants.H"
#include "BOOM/GffReader.H"
#include "BOOM/Regex.H"
#include "EmissionLoader.H"
#include "SequenceSet.H"
using namespace std;

class Application {
  Map<String,Vector<GffFeature*> > intervalsBySubstrate;
  int trackID;
  Regex seqfile;
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
  : seqfile("([^\/]+)\\.fastb")
{
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"l");
  if(cmd.numArgs()!=4)
    throw string(
"\n\
motif-enrichment <fastb-dir> <schema> <trackname> <comma,separated,motif,list>\n\
    where: -l = report logarithmic value\n\
\n\
"
);
  String fastbDir=cmd.arg(0);
  String schemaFile=cmd.arg(1);
  String trackName=cmd.arg(2);
  String motifList=cmd.arg(3);
  BOOM::Vector<BOOM::String> &motifs=*motifList.getFields(",");
  const int numMotifs=motifs.size();
  bool wantLog=cmd.option('l');

  // Load sequences
  Schema schema(schemaFile);
  const int trackID=schema.lookupDiscreteID(trackName);
  EmissionLoader loader(schema);
  Vector<String> files;
  File::getFileList(fastbDir,files);
  int numFiles=files.size();
  int totalHits=0, totalDNA=0;
  for(int k=0 ; k<numFiles ; ++k) {
    const String &filename=files[k];
    if(!seqfile.search(filename)) continue;
    String substrate=seqfile[1];
    String path=fastbDir+"/"+filename;
    EmissionSequence *seq=loader.load(path);
    const int seqLen=seq->length();
    totalDNA+=seqLen;
    Sequence &dnaSeq=*seq->getDiscreteSeq(trackID);
    String dna=dnaSeq(schema.getAlphabet(trackID));
    for(int i=0 ; i<numMotifs ; ++i) {
      String &motif=motifs[i];
      int motifLen=motif.size();
      int lastPos=seqLen-motifLen;
      for(int pos=0 ; pos<=lastPos ; ++pos) 
	if(dna.occursAt(motif,pos)) ++totalHits;
    }
    delete &dnaSeq;
    delete seq;
  }
  float density=totalHits/float(totalDNA);
  if(wantLog) density=log(totalHits)-log(totalDNA);
  cout.precision(30);
  cout<<density<<endl;

  return 0;
}



