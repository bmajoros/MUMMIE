/****************************************************************
 fastb-slice.C
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Interval.H"
#include "EmissionLoader.H"
using namespace std;
using namespace BOOM;

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
};


int main(int argc,char *argv[])
{
  try {
    Application app;
    return app.main(argc,argv);
  }
  catch(const char *p) { cerr << p << endl; }
  catch(const string &msg) { cerr << msg.c_str() << endl; }
  catch(const exception &e)
    {cerr << "STL exception caught in main:\n" << e.what() << endl;}
  catch(...) { cerr << "Unknown exception caught in main" << endl; }
  return -1;
}



Application::Application()
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=5)
    throw String("\n\
fastb-slice <in.schema> <in.fastb> <begin> <end> <out.fastb>\n\
");
  const String schemaFile=cmd.arg(0);
  const String infile=cmd.arg(1);
  const int begin=cmd.arg(2).asInt();
  const int end=cmd.arg(3).asInt();
  const String outFile=cmd.arg(4);

  // Load the sequence and schema
  Schema schema(schemaFile);
  EmissionLoader loader(schema);
  EmissionSequence &S=*loader.load(infile);

  // Slice the fastb
  EmissionSequence *sub=S.getSubsequence(begin,end-begin);
  sub->save(outFile);
  delete sub;

  // Clean up
  delete &S;
  return 0;
}




