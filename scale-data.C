/**************************************************************
 scale-data.C 
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#include <string>
#include <iostream>
#include <fstream>
#include "BOOM/Random.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Vector.H"
#include "BOOM/SparseGraph.H"
#include "BOOM/Constants.H"
#include "EmissionLoader.H"
#include "SequenceSet.H"
#include "Schema.H"
using namespace std;

class Application
{
  SequenceSet trainingSet;
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
  CommandLine cmd(argc,argv,"n:");
  if(cmd.numArgs()!=3) 
    throw "\n\
scale-data [options] <input-dir> <transforms.out> <schema.txt>\n\
   where: -n V = add gaussian noise (after scaling) with variance V\n\
\n\
";
  String inDir=cmd.arg(0);
  String transformFile=cmd.arg(1);
  String schemaFile=cmd.arg(2);
  bool wantNoise=cmd.option('n');
  double noiseVar=wantNoise ? cmd.optParm('n').asDouble() : 0;

  // Misc. initialization
  cout.precision(10);
  cerr.precision(10);
  if(cmd.option('s')) SeedRandomizer(cmd.optParm('s').asUnsigned());
  else randomize();
  BOOM::catchFloatOverflow();

  // Load the training set
  Schema schema(schemaFile);
  cout<<"loading training data"<<endl;
  trainingSet.load(inDir,schema);

  // Rescale the training data
  cout<<"rescaling data..."<<endl;
  ofstream os(transformFile.c_str());
  trainingSet.normalize();
  trainingSet.saveTransforms(os);
  os.close();

  // Add noise
  if(wantNoise) trainingSet.addNoise(0,noiseVar);

  // Save the results
  trainingSet.save();

  return 0;
}


