/**************************************************************
 threshold.C 
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#include <string>
#include <iostream>
#include <fstream>
#include "BOOM/Random.H"
#include "BOOM/CommandLine.H"
#include "BOOM/VectorSorter.H"
#include "BOOM/Array1D.H"
#include "BOOM/Array2D.H"
#include "BOOM/SparseGraph.H"
#include "BOOM/Constants.H"
#include "BOOM/Exceptions.H"
#include "EmissionLoader.H"
#include "SequenceSet.H"
#include "Schema.H"
using namespace std;

class Application {
  Schema schema;
  int D;
  double threshold;
  SequenceSet trainingSet;
  void applyThreshold();
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
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=3)
    throw "\n\
threshold <in.schema> <input-dir> <threshold>\n\
\n\
";
  String schemaFile=cmd.arg(0);
  String inDir=cmd.arg(1);
  threshold=cmd.arg(2).asDouble();

  // Misc. initialization
  cout.precision(10);
  cerr.precision(10);
  BOOM::catchFloatOverflow();

  // Load the training set
  schema=Schema(schemaFile);
  D=schema.getNumContinuous();
  trainingSet.load(inDir,schema);

  // Apply thresholding to the data
  applyThreshold();

  // Save the results
  trainingSet.save();

  return 0;
}



void Application::applyThreshold()
{
  const int numSeqs=trainingSet.size();
  for(int i=0 ; i<numSeqs ; ++i) {
    EmissionSequence &S=*trainingSet[i];
    const int L=S.length();
    for(int pos=0 ; pos<L ; ++pos) {
      Emission &e=S[pos];
      GSL::Vector &V=e.getContinuous();
      for(int d=0 ; d<D ; ++d) V[d]=V[d]>=threshold ? threshold : 0;
    }
  }
}




