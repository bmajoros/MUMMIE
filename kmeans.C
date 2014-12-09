/**************************************************************
 kmeans.C
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
#include "BOOM/Constants.H"
#include "BOOM/Random.H"
#include "EmissionLoader.H"
#include "Kmeans.H"
#include "Schema.H"
using namespace std;

class Application {
  BOOM::Vector<GSL::Vector> data;
  Schema *schema;
  void loadData(String dirName,int maxCases);
  void transpose();
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
  CommandLine cmd(argc,argv,"n:i:R:r:s:tmc:");
  if(cmd.numArgs()!=3)
    throw "\n\
kmeans [options] <training-dir> <K> <*.schema>\n\
   where:  -n NUM = use at most NUM training cases\n\
           -i NUM = perform at most NUM iterations\n\
           -R epsilon = stop when change in RSS<epsilon\n\
           -r NUM = use NUM random restarts; keep best RSS solution\n\
           -s NUM = use NUM iterations of random swapping upon convergence\n\
           -t = tranpose the data matrix before clustering\n\
           -m = print membership table\n\
           -c file = write centroids into file\n\
";
  String infile=cmd.arg(0);
  const int K=cmd.arg(1).asInt();
  const String schemaFile=cmd.arg(2);
  const int maxCases=cmd.option('n') ? cmd.optParm('n').asInt() : -1;
  const int maxIter=cmd.option('i') ? cmd.optParm('i').asInt() : -1;
  const int numRestarts=cmd.option('r') ? cmd.optParm('r').asInt() : 0;
  const int numSwaps=cmd.option('s') ? cmd.optParm('s').asInt() : 0;
  const double minChange=cmd.option('R') ? cmd.optParm('R').asDouble() : -1;
  const bool shouldTranspose=cmd.option('t');
  const bool wantMembership=cmd.option('m');
  const String outfile=cmd.option('c') ? cmd.optParm('c') : "";
  if(numRestarts>0 && numSwaps>0) throw "-s and -r cannot be used together";

  // Misc. initialization
  randomize();
  BOOM::catchFloatOverflow();
  schema=new Schema(schemaFile);

  // Load the data
  loadData(infile,maxCases);
  if(shouldTranspose) transpose();

  // Perform clustering
  KmeansEuclidean cmp;
  Kmeans kmeans(K,cmp);
  double RSS;
  if(numSwaps>0) RSS=kmeans.cluster_swap(data,maxIter,minChange,numSwaps);
  else RSS=kmeans.cluster(data,maxIter,minChange,numRestarts);
  const Array1D<GSL::Vector> &centroids=kmeans.getCentroids();
  cout<<"best RSS="<<RSS<<endl;
  if(!outfile.isEmpty()) {
    ofstream os(outfile.c_str());
    os.precision(30);
    os<<K<<endl;
    for(int i=0 ; i<K ; ++i) os<<centroids[i]<<endl;
  }
  if(wantMembership) {
    cout<<"cluster sizes: ";
    const IntArray1D &M=kmeans.getMembership();
    IntArray1D counts(K); counts.setAllTo(0);
    const int n=M.size();
    for(int i=0 ; i<n ; ++i) ++counts[M[i]];
    for(int i=0 ; i<K ; ++i) cout<<counts[i]<<" ";
    cout<<endl;
  }

  return 0;
}


void Application::loadData(String dir,int maxCases)
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
  EmissionLoader loader(*schema);
  for(int i=0 ; i<numFiles ; ++i) {
    EmissionSequence &s=*loader.load(dir+"/"+files[i]);
    const int L=s.length();
    for(int j=0 ; j<L ; ++j) {
      data.push_back(s[j].getContinuous());
      if(maxCases>0 && data.size()>=maxCases) break;
    }
    delete &s;
    if(maxCases>0 && data.size()>=maxCases) break;
  }
}



void Application::transpose()
{
  BOOM::Vector<GSL::Vector> T;
  int dim=data[0].getDim();
  int N=data.size();
  for(int i=0 ; i<dim ; ++i) {
    GSL::Vector rec(N);
    for(int j=0 ; j<N ; ++j)
      rec[j]=data[j][i];
    T.push_back(rec);
  }
  data=T;
}


