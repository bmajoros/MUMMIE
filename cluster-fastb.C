/****************************************************************
 cluster-fastb.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/Map.H"
#include "BOOM/Regex.H"
#include "BOOM/CommandLine.H"
#include "BOOM/SummaryStats.H"
#include "BOOM/Correlation.H"
#include "BOOM/GSL/Matrix.H"
#include "BOOM/GSL/Vector.H"
#include "EmissionLoader.H"
using namespace std;
using namespace BOOM;


class Application {
public:
  Application();
  int main(int argc,char *argv[]);
protected:
  Regex regex;
  int D;
  GSL::Matrix cov, inverseCov;
  GSL::Vector means;
  BOOM::Vector<EmissionSequence*> sequences;
  void loadSequences(String dirName);
  double mahalanobis(GSL::Vector &,GSL::Vector &);
  void Application::computeCov();
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
  : regex("(.*)...-Dmel")
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=2)
    throw String("cluster-fastb <fastb-directory> <cutoff>");
  const String dir=cmd.arg(0);
  const double cutoff=cmd.arg(1).asDouble();

  Map<String,int> lengths;
  ifstream is("lengths.txt");
  while(!is.eof()) {
    String line;
    line.getline(is);
    BOOM::Vector<String> &fields=*line.getFields();
    if(&fields!=NULL && fields.size()>1) {
      lengths[fields[0]]=fields[1].asInt();
      //cout<<fields[0]<<" -> "<<fields[1].asInt()<<endl;
    }
    delete &fields;
  }

  loadSequences(dir);
  D=(*sequences[0])[0].getContinuous().getDim();
  /*
  computeCov();
  cov.invert(inverseCov);
  */

  int numSeqs=sequences.size();
  for(int j=0 ; j<numSeqs ; ++j) {
    EmissionSequence &seq=*sequences[j];
    int L=seq.length();
    //cout<<"L="<<L<<endl;
    const Schema &schema=seq.getSchema();
    String outfile=String("for-iulian/")+(j+1)+".txt";
    ofstream os(outfile.c_str());
    for(int d=0 ; d<D ; ++d) {
      const String &name=schema.getContinuousName(d);
      if(regex.match(name)) name=regex[1];
      const int len=lengths[name];
      if(len==0 || len>20) {cout<<len<<" "<<name<<endl;INTERNAL_ERROR;}
      for(int i=1 ; i<L ; ++i) {
	double x=seq[i].getContinuous()[d];
	for(int j=i-1 ; i-j<len && j>=0 ; --j) {
	  double t=seq[j].getContinuous()[d];
	  if(t>=x) x=t; 
	}
	if(x<cutoff) x=0.0;
	os<<x<<"\t";
      }
      os<<endl;
    }
  }
  /*
  BOOM::Vector<EmissionSequence*>::iterator cur=sequences.begin(), 
    end=sequences.end();
  for(; cur!=end ; ++cur) {
    EmissionSequence &seq=**cur;
    int L=seq.length();
    for(int i=0 ; i<L ; ++i) {
      const Emission &e=seq[i];
      const GSL::Vector &v=e.getContinuous();
      for(int d=0 ; d<D ; ++d) {
	cout<<v[d]<<"\t";
      }
      cout<<endl;
    }
  }
   */
  return 0;
}



double Application::mahalanobis(GSL::Vector &v1,
				GSL::Vector &v2)
{
  int D=v1.getDim();
  GSL::Matrix left(1,D), right(D,1);
  for(int i=0 ; i<D ; ++i) {
    double diff=v1[i]-v2[i];
    left(0,i)=right(i,0)=diff;
  }
  GSL::Matrix temp1(1,D), temp2(1,1);
  left.times(inverseCov,temp1);
  temp1.times(right,temp2);
  double mah=temp2(0,0);
  return mah;
}



void Application::computeCov()
{
  means.resize(D);
  cov.resize(D,D);
  const int numSeqs=sequences.size();

  // Compute means and variances
  for(int d=0 ; d<D ; ++d) {
    DoubleVector V;
    for(int i=0 ; i<numSeqs ; ++i) {
      const EmissionSequence &seq=*sequences[i];
      const int L=seq.length();
      for(int pos=0 ; pos<L ; ++pos) {
	const Emission &s=seq[pos];
	const GSL::Vector &v=s.getContinuous();
	V.push_back(v[d]);
      }
    }
    SummaryStats stats(V);
    means[d]=stats.getMean();
    cov(d,d)=stats.getVar();
  }
  
  // Compute covariances
  for(int d1=0 ; d1<D ; ++d1) {
    double sd1=sqrt(cov(d1,d1));
    for(int d2=0 ; d2<D ; ++d2) {
      if(d2==d1) continue;
      double sd2=sqrt(cov(d2,d2));
      DoubleVector V1, V2;
      for(int i=0 ; i<numSeqs ; ++i) {
	const EmissionSequence &seq=*sequences[i];
	const int L=seq.length();
	for(int pos=0 ; pos<L ; ++pos) {
	  const Emission &s=seq[pos];
	  const GSL::Vector &v=s.getContinuous();
	  V1.push_back(v[d1]);
	  V2.push_back(v[d2]);
	}
      }
      Correlation cor(V1,V2);
      cov(d1,d2)=cor.getR()*sd1*sd2;
    }
  }
}



void Application::loadSequences(String dir)
{
  BOOM::Vector<String> files, tmp;
  File::getFileList(dir,files);
  int numFiles=files.size();
  for(int i=0 ; i<numFiles ; ++i) {
    String filename=files[i];
    int L=filename.length();
    if(L<7) continue;
    String extension=filename.substr(L-5);
    extension.toupper();
    if(extension=="FASTB") tmp.push_back(filename);
  }
  files=tmp;
  numFiles=files.size();
  if(numFiles==0) throw "No files with extension \".fastb\" found";
  EmissionLoader loader;
  for(int i=0 ; i<numFiles ; ++i) {
    cout<<"loading "<<files[i]<<endl;
    sequences.push_back(loader.load(dir+"/"+files[i]));
  }
}


