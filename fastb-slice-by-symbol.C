/****************************************************************
 fastb-slice-by-symbol.C
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
private:
  void findIntervals(const EmissionSequence &,
		     int trackID,
		     int symbol,
		     Vector<Interval> &);
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
  CommandLine cmd(argc,argv,"d");
  if(cmd.numArgs()!=6)
    throw String("\n\
fastb-slice-by-symbol <in.schema> <in.fastb> <out-dir> <out-filestem> <track> <symbol>\n\
   -d : drop the discrete track when writing output files\n\
");
  const String schemaFile=cmd.arg(0);
  const String infile=cmd.arg(1);
  const String outDir=cmd.arg(2);
  const String filestem=cmd.arg(3);
  const String trackName=cmd.arg(4);
  const String symbolStr=cmd.arg(5);
  const bool wantDrop=cmd.option('d');
  if(symbolStr.length()!=1) throw "symbol must be a single character";
  if(outDir.length()>0 && outDir.lastChar()=='/') outDir.chop();

  // Load the sequence and schema
  Schema schema(schemaFile);
  EmissionLoader loader(schema);
  EmissionSequence &S=*loader.load(infile);
  int trackID=schema.lookupDiscreteID(trackName);
  Alphabet &alpha=schema.getAlphabet(trackID);
  int symbol=alpha.lookup(symbolStr[0]);

  // Find intervals marked by the symbol and emit them
  Vector<Interval> intervals;
  findIntervals(S,trackID,symbol,intervals);
  int nextFile=1;
  for(Vector<Interval>::const_iterator cur=intervals.begin(), 
	end=intervals.end() ; cur!=end ; ++cur) {
    const Interval I=*cur;
    EmissionSequence *sub=S.getSubsequence(I.getBegin(),I.length());
    if(wantDrop) sub->dropDiscreteTrack(trackID);
    String filename=outDir+"/"+filestem+nextFile+".fastb";
    sub->save(filename);
    delete sub;
    ++nextFile;
  }

  // Clean up
  delete &S;
  return 0;
}



void Application::findIntervals(const EmissionSequence &S,
				int trackID,
				int target,
				Vector<Interval> &intervals)
{
  const int L=S.length();
  int begin=-1;
  for(int i=0 ; i<L ; ++i) {
    NmerSymbol s=S[i].getDiscrete(trackID);
    if(s==target) { if(begin<0) begin=i; }
    else if(begin>=0) {
      intervals.push_back(Interval(begin,i));
      begin=-1;
    }
  }
  if(begin>=0) intervals.push_back(Interval(begin,L));
}



