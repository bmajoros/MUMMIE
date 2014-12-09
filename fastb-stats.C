/**************************************************************
 fastb-stats.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#include <string>
#include <iostream>
#include <fstream>
#include "BOOM/CommandLine.H"
#include "BOOM/Vector.H"
#include "BOOM/Array2D.H"
#include "BOOM/Constants.H"
#include "BOOM/SummaryStats.H"
#include "BOOM/Correlation.H"
#include "EmissionLoader.H"
using namespace std;

class Application {
  void update(EmissionSequence &,Schema &,Array1D<double> &sums,
	      Array1D<double> &squares,Array2D<double> &products,int &N);
  void report(EmissionSequence &,Schema &);
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
  CommandLine cmd(argc,argv,"ds");
  if(cmd.numArgs()!=2) 
    throw "\n\
fastb-stats [options] <*.fastb> <*.schema>\n\
  where:  -d = filename is actually a directory\n\
          -s = give overall summary for whole directory (not each file)\n\
\n\
";
  String infile=cmd.arg(0);
  String schemaFile=cmd.arg(1);
  bool isDir=cmd.option('d');
  //cout.precision(30);
  //cerr.precision(30);

  // Misc. initialization
  BOOM::catchFloatOverflow();

  // Process the fastb file(s)
  Schema schema(schemaFile);
  EmissionLoader loader(schema);
  if(cmd.option('s')) {
    int numC=schema.getNumContinuous();
    Array1D<double> sums(numC); sums.setAllTo(0.0);
    Array1D<double> squares(numC); squares.setAllTo(0.0);
    Array2D<double> products(numC,numC); products.setAllTo(0.0);
    int N=0;
    const String dir=infile;
    BOOM::Vector<String> files;
    File::getFileList(dir,files);
    int numFiles=files.size();
    for(int i=0 ; i<numFiles ; ++i) {
      String filename=files[i];
      int L=filename.length();
      if(L<7) continue;
      String extension=filename.substr(L-5);
      extension.toupper();
      if(extension!="FASTB") continue;
      EmissionSequence *S=loader.load(dir+"/"+filename);
      update(*S,schema,sums,squares,products,N);
      delete S;
    }
    for(int i=0 ; i<numC ; ++i) {
      String name=schema.getContinuousName(i);
      const double sumX=sums[i], sumXX=squares[i];
      double mean=sumX/N;
      double stddev=sqrt((sumXX-sumX*sumX/N)/(N-1.0));
      cout<<name<<"\t"<<mean<<"\t"<<stddev<<" N="<<N<<endl;
    }
    for(int i=0 ; i<numC ; ++i) {
      String name=schema.getContinuousName(i);
      for(int j=i+1 ; j<numC ; ++j) 
	cout<<name<<":"<<schema.getContinuousName(j)<<" r="<<
	  Correlation::computeR(sums[i],products[i][i],sums[j],products[j][j],
				products[i][j],N)<<endl;
    }
    return 0;
  }
  if(isDir) {
    const String dir=infile;
    BOOM::Vector<String> files;
    File::getFileList(dir,files);
    int numFiles=files.size();
    for(int i=0 ; i<numFiles ; ++i) {
      String filename=files[i];
      int L=filename.length();
      if(L<7) continue;
      String extension=filename.substr(L-5);
      extension.toupper();
      if(extension!="FASTB") continue;
      EmissionSequence *S=loader.load(dir+"/"+filename);
      report(*S,schema);
      delete S;
    }
  }
  else {
    EmissionSequence *S=loader.load(infile);
    report(*S,schema);
    delete S;
  }
  
  return 0;
}


void Application::report(EmissionSequence &S,Schema &schema)
{
  int len=S.length();
  cout<<"length "<<len<<endl;
  int numD=schema.getNumDiscrete();
  int numC=schema.getNumContinuous();
  for(int i=0 ; i<numD ; ++i) {
    const String &name=schema.getDiscreteName(i);
    cout<<name<<" ";
    Sequence *seq=S.getDiscreteSeq(i);
    Alphabet &alpha=schema.getAlphabet(i);
    int numAlpha=alpha.getNumElements();
    for(Symbol s=0 ; s<numAlpha ; ++s) {
      int n=seq->countOccurrences(s);
      float percent=int(n/float(len)*1000+5/9.0)/10.0;
      cout<<alpha.lookup(s)<<"="<<percent<<"% ";
    }
    int n=seq->countOccurrences(-1);
    float percent=int(n/float(len)*1000+5/9.0)/10.0;
    cout<<"?="<<percent<<"% ";
    cout<<endl;
    delete seq;
  }
  for(int i=0 ; i<numC ; ++i) {
    const String &name=schema.getContinuousName(i);
    BOOM::Vector<double> V;
    for(int pos=0 ; pos<len ; ++pos) V.push_back(S[pos][i]);
    SummaryStats stats(V);
    cout<<name<<" "<<stats.getMean()<<" +/- "<<stats.getStdDev()
	<<" ("<<stats.getMin()<<" - "<<stats.getMax()<<")"<<endl;
  }
}



void Application::update(EmissionSequence &S,Schema &schema,
			 Array1D<double> &sums,
			 Array1D<double> &squares,
			 Array2D<double> &products,int &N)
{
  int len=S.length();
  int numC=schema.getNumContinuous();
  for(int pos=0 ; pos<len ; ++pos)
    for(int i=0 ; i<numC ; ++i) {
      Emission &E=S[pos];
      const double X=E[i];
      sums[i]+=X;
      squares[i]+=X*X;
      for(int j=0 ; j<numC ; ++j) 
	products[i][j]+=X*E[j];
    }
  N+=len;
}


