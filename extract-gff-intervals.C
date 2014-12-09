/*
 extract-gff-intervals.C
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
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=4)
    throw string(
"\n\
extract-gff-intervals <in.gff> <fastb-dir> <schema> <out-dir>\n\
\n\
"
);
  String gffFile=cmd.arg(0);
  String fastbDir=cmd.arg(1);
  String schemaFile=cmd.arg(2);
  String outdir=cmd.arg(3);

  // Load GFF
  GffReader reader(gffFile);
  Vector<GffFeature*> &features=*reader.loadFeatures();
  int numFeatures=features.size();
  for(int i=0 ; i<numFeatures ; ++i) {
    GffFeature &feature=*features[i];
    const String &substrate=feature.getSubstrate();
    intervalsBySubstrate[substrate].push_back(&feature);
  }

  // Load sequences
  Schema schema(schemaFile);
  EmissionLoader loader(schema);
  Vector<String> files;
  File::getFileList(fastbDir,files);
  int numFiles=files.size();
  for(int k=0 ; k<numFiles ; ++k) {
    const String &filename=files[k];
    if(!seqfile.search(filename)) continue;
      //throw String("Can't parse filename: ")+filename;
    String substrate=seqfile[1];
    String path=fastbDir+"/"+filename;
    if(!File::exists(path)) continue;
    EmissionSequence *seq=loader.load(path);
    const int seqLen=seq->length();
    if(intervalsBySubstrate.isDefined(substrate)) {
      Vector<GffFeature*> &intervals=intervalsBySubstrate[substrate];
      int numIntervals=intervals.size();
      for(int i=0 ; i<numIntervals ; ++i) {
	GffFeature *feature=intervals[i];
	int begin=feature->getBegin(), end=feature->getEnd();
	int length=end-begin;
	if(length<1) {
	  cout<<*feature<<endl;
	  throw String("bad length: begin=")+begin+" end="+end+" L="+seqLen;
	}
	EmissionSequence *sub=seq->getSubsequence(begin,length);
	String outpath=outdir+"/"+substrate+"-interval"+(i+1)+".fastb";
	sub->save(outpath);
	delete sub;
      }
    }
    delete seq;
  }

  return 0;
}



