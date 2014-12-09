/****************************************************************
 subset-fastb.C
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
using namespace std;
using namespace BOOM;


class Application {
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
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=3)
    throw String("subset-fastb <in.fastb> <track,track,...> <out.fastb>");
  const String infile=cmd.arg(0);
  const String trackListStr=cmd.arg(1);
  const String outfile=cmd.arg(2);

  BOOM::Vector<String> &fields=*trackListStr.getFields(",");
  Set<String> keep;
  int n=fields.size();
  for(int i=0 ; i<n ; ++i) keep.insert(fields[i]);

  ofstream os(outfile.c_str());
  FastbReader reader(infile);
  while(!reader.eof()) {
    FastbSequence *seq=reader.nextSequence();
    if(!seq) break;
    if(keep.isMember(seq->getID())) os<<*seq<<endl;
    delete seq;
  }	

  return 0;
}






