/**************************************************************
fastb-add-wig-track.C
W.H. Majoros
bmajoros@duke.edu

This is OPEN SOURCE SOFTWARE governed by the ARTISTIC LICENSE.
***************************************************************/
#include <string>
#include <iostream>
#include <fstream>
#include "BOOM/CommandLine.H"
#include "BOOM/Vector.H"
#include "BOOM/Constants.H"
#include "BOOM/SummaryStats.H"
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
  if(cmd.numArgs()!=4)
    throw "\n\
fastb-add-wig-track <schema> <indir> <wigdir> <outdir>\n\
\n\
";
  String schemqaFile=cmd.arg(0);
  String inDir=cmd.arg(1);
  String wigDir=cmd.arg(2);
  String outDir=cmd.arg(3);
  //cout.precision(30);
  //cerr.precision(30);

  // Misc. initialization
  BOOM::catchFloatOverflow();

  // Process the fastb file(s)
  Schema schema(schemaFile);
  EmissionLoader loader(schema);
  BOOM::Vector<String> files;
  File::getFileList(inDir,files);
  int numFiles=files.size();
  for(int i=0 ; i<numFiles ; ++i) {
    String filename=files[i];
    int L=filename.length();
    if(L<7) continue;
    String extension=filename.substr(L-5);
    extension.toupper();
    if(extension!="FASTB") continue;
    EmissionSequence *S=loader.load(dir+"/"+filename);


    delete S;
  }
  
  return 0;
}



