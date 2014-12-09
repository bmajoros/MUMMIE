/**************************************************************
 discretize.C 
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
  Schema schema, newSchema;
  int D, numDiscrete;
  int numBins;
  SequenceSet trainingSet;
  Array2D<double> boundaries; // [track][boundary]
  void analyze();
  void discretize();
  NmerSymbol discretize(double value,int dim);
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
  if(cmd.numArgs()!=5)
    throw "\n\
discretize <in.schema> <in/out-dir> <num-bins> <boundaries.out> <out.schema>\n\
\n\
";
  String schemaFile=cmd.arg(0);
  String inDir=cmd.arg(1);
  numBins=cmd.arg(2).asInt();
  String transformFile=cmd.arg(3);
  String schemaOutFile=cmd.arg(4);

  // Misc. initialization
  cout.precision(10);
  cerr.precision(10);
  BOOM::catchFloatOverflow();

  // Load the training set
  schema=Schema(schemaFile);
  D=schema.getNumContinuous();
  numDiscrete=schema.getNumDiscrete();
  trainingSet.load(inDir,schema);

  // Create new schema
  const char firstLetter='A';
  newSchema=Schema(D+numDiscrete);
  for(int d=0 ; d<D ; ++d) {
    Alphabet &A=newSchema.getAlphabet(d);
    newSchema.getDiscreteName(d)=schema.getContinuousName(d);
    for(int s=0 ; s<numBins ; ++s) A.add(firstLetter+s);
  }
  for(int d=0 ; d<numDiscrete ; ++d) {
    newSchema.getDiscreteName(d+D)=schema.getDiscreteName(d);
    newSchema.getAlphabet(d+D)=schema.getAlphabet(d);
  }	

  // Analyze the data
  analyze();

  // Write transform to file
  ofstream os(transformFile.c_str()), os2(schemaOutFile.c_str());
  for(int d=0 ; d<D ; ++d) {
    for(int b=0 ; b<numBins-1 ; ++b)
      os<<boundaries[d][b]<<"\t";
    os<<endl;
  }
  os.close();
  newSchema.save(os2);
  os2.close();

  // Discretize
  discretize();

  // Save the results
  trainingSet.save();

  return 0;
}



void Application::analyze()
{
  Array1D< Vector<double> > data(D);
  boundaries.resize(D,numBins-1);
  const int numSeqs=trainingSet.size();
  for(int i=0 ; i<numSeqs ; ++i) {
    EmissionSequence &S=*trainingSet[i];
    const int L=S.length();
    for(int pos=0 ; pos<L ; ++pos) {
      Emission &e=S[pos];
      for(int d=0 ; d<D ; ++d) data[d].push_back(e[d]);
    }
  }
  DirectComparator<double> cmp;
  for(int d=0 ; d<D ; ++d) {
    Vector<double> &V=data[d];
    const int L=V.size();
    VectorSorter<double> sorter(V,cmp);
    sorter.sortAscendInPlace();
    Vector<double> &v=*V.uniq();
    const int l=v.size();
    for(int b=0 ; b<numBins-1 ; ++b)
      boundaries[d][b]=v[(b+1)*l/double(numBins)];
    delete &v;
  }
}



void Application::discretize()
{
  const int numSeqs=trainingSet.size();
  for(int i=0 ; i<numSeqs ; ++i) {
    EmissionSequence &S=*trainingSet[i], T(newSchema);
    const int L=S.length();
    for(int pos=0 ; pos<L ; ++pos) {
      Emission &e=S[pos], f;
      f.resize(D+numDiscrete,0);
      GSL::Vector &V=e.getContinuous();
      for(int d=0 ; d<D ; ++d)
	f.getDiscrete(d)=discretize(V[d],d);
      for(int d=0 ; d<numDiscrete ; ++d)
	f.getDiscrete(D+d)=e.getDiscrete(d);
      T.append(f);
    }
    S.copySeqAndSchema(T);
  }
}



NmerSymbol Application::discretize(double value,int dim)
{
  Array2D<double>::RowIn2DArray<double> row=boundaries[dim];
  for(int b=0 ; b<numBins-1 ; ++b)
    if(value<row[b]) return b;
  return numBins-1;
}



