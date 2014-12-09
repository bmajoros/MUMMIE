/*
 viterbi-trainer.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
*/
#include <math.h>
#include <iostream>
#include "BOOM/CommandLine.H"
#include "BOOM/Vector.H"
#include "BOOM/Array2D.H"
#include "BOOM/Regex.H"
#include "BOOM/File.H"
#include "BOOM/Correlation.H"
#include "BOOM/SummaryStats.H"
#include "ComponentViterbi.H"
#include "EmissionSequence.H"
#include "EmissionLoader.H"
using namespace std;


struct EmissionMatrix {
  EmissionMatrix(int d1,int d2,Schema &schema) :M(d1,d2), d1(d1), d2(d2) {
    for(int i=0;i<d1;++i) for(int j=0;j<d2;++j)
      M[i][j]=new EmissionSequence(schema);
  }
  virtual ~EmissionMatrix() {
    for(int i=0;i<d1;++i) for(int j=0;j<d2;++j) delete M[i][j];
  }
  EmissionSequence &operator()(int i,int j) {return *M[i][j];}
protected:
  int d1, d2;
  Array2D<EmissionSequence*> M;
};



class Application {
  bool wantPoolStates;
  BOOM::Vector<EmissionSequence*> trainingSet;
  void loadTrainingSet(String dirName);
  void reestimate(EmissionSequence &,MultiGauss &);
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
  CommandLine cmd(argc,argv,"s");
  if(cmd.numArgs()!=4)
    throw string(
"viterbi-trainer [-s] <num-iterations> <template.hmm> <training-dir>\n\
                      <outfile.hmm>\n\
   where: -s = pool emission distributions across states\n\
\n\
");
  int numIterations=cmd.arg(0).asInt();
  String inHmmFile=cmd.arg(1);
  String dir=cmd.arg(2);
  String outfile=cmd.arg(3);
  wantPoolStates=cmd.option('s');
  
  // Load HMM and perform some initialization
  HMM hmm(inHmmFile);
  ComponentViterbi viterbi(hmm);
  int numStates=hmm.countStates();
  int numComponents=hmm.numMixtureComponents();
  Schema &schema=hmm.getSchema();
  
  // Load training set
  loadTrainingSet(dir);
  int numTrain=trainingSet.size();
  
  // Iterate through all training sequences
  BOOM::Vector<StateComponentPair> prevPath;
  bool done=false;
  for(int k=0 ; k<numIterations && !done ; ++k) {
    Array2D<float> transCounts(numStates,numStates);
    transCounts.setAllTo(0.0);
    EmissionMatrix emissions(numStates,numComponents,schema);
    for(int i=0 ; i<numTrain ; ++i) {
      cout<<"iteration "<<k<<endl;
      EmissionSequence &seq=*trainingSet[i];
      double pathScore;
      BOOM::Vector<StateComponentPair> &path=*viterbi.getPath(seq,pathScore);
      if(path==prevPath) {done=true; break;}
      prevPath=path;
      path.push_front(StateComponentPair(0,0));
      path.push_back(StateComponentPair(0,0));
      //cout<<"path="<<path<<endl;
      int len=path.size();
      for(int j=0 ; j<len-1 ; ++j){
	int fromQ=path[j].first, toQ=path[j+1].first;
	cout<<fromQ<<" ";
	++transCounts[fromQ][toQ];
      }
      cout<<endl;
      for(int j=1 ; j<len-1 ; ++j) {
	const StateComponentPair &p=path[j];
	emissions(p.first,p.second).append(seq[j-1]);
      }
      delete &path;
    }
    if(done) break;

    // Normalize transition counts
    for(int i=0 ; i<numStates ; ++i) {
      double sum=0;
      for(int j=0 ; j<numStates ; ++j) {
	sum+=transCounts[i][j];
      }
      for(int j=0 ; j<numStates ; ++j) {
	double p=transCounts[i][j]/sum;
	hmm.setTransitionProb(i,j,p);
      }
    }

    // Reestimate emissions
    for(int q=1 ; q<numStates ; ++q) {
      GaussianMixture &mix=hmm.getEmissionDistr(q);
      int numComponents=mix.getNumComponents();
      Array1D<int> coefCounts(numComponents);
      coefCounts.setAllTo(0);
      for(int c=0 ; c<numComponents ; ++c) 
	coefCounts[c]+=emissions(q,c).length();
      double sum=0;
      for(int c=0 ; c<numComponents ; ++c) sum+=coefCounts[c];
      for(int c=0 ; c<numComponents ; ++c)
	mix.setCoef(c,coefCounts[c]/sum);
      if(!wantPoolStates)
	for(int c=0 ; c<numComponents ; ++c) {
	  MultiGauss newDistr;
	  EmissionSequence &seq=emissions(q,c);
	  if(seq.length()==0) throw "Component has no data";
	  reestimate(seq,newDistr);
	  mix.setDistr(c,newDistr);
	}
    }
    if(wantPoolStates)
      cout<<"POOLING"<<endl;
      for(int c=0 ; c<numComponents ; ++c) {
	MultiGauss newDistr;
	EmissionSequence pooled(schema);
	for(STATE q=0 ; q<numStates ; ++q) {
	  EmissionSequence &seq=emissions(q,c);
	  pooled.append(seq);
	}
	if(pooled.length()==0) throw "Component has no data";
	reestimate(pooled,newDistr);
	for(int q=1 ; q<numStates ; ++q) {
	  GaussianMixture &mix=hmm.getEmissionDistr(q);
	  mix.setDistr(c,newDistr);
	}
      }
    cout<<hmm<<endl;
  }
  hmm.save(outfile);  
  return 0;
}



void Application::reestimate(EmissionSequence &seq,MultiGauss &dist)
{
  int dim=seq.getSchema().getNumContinuous();
  GSL::Vector means(dim);
  GSL::Matrix covM(dim,dim);
  int L=seq.length();
  for(int d=0 ; d<dim ; ++d) {
    DoubleVector V;
    for(int pos=0 ; pos<L ; ++pos) {
      const Emission &s=seq[pos];
      const GSL::Vector &v=s.getContinuous();
      V.push_back(v[d]);
    }
    SummaryStats stats(V);
    means[d]=stats.getMean();
    covM(d,d)=stats.getVar();
  }
  for(int d1=0 ; d1<dim ; ++d1) {
    double sd1=sqrt(covM(d1,d1));
    for(int d2=0 ; d2<dim ; ++d2) {
      if(d2==d1) continue;
      double sd2=sqrt(covM(d2,d2));
      DoubleVector V1, V2;
      for(int pos=0 ; pos<L ; ++pos) {
	const Emission &s=seq[pos];
	const GSL::Vector &v=s.getContinuous();
	V1.push_back(v[d1]);
	V2.push_back(v[d2]);
      }
      Correlation cor(V1,V2);
      covM(d1,d2)=cor.getR()*sd1*sd2;
    }
  }
  cout<<"REESTIMATING FROM: "<<means<<"\n"<<covM<<endl;
  dist=MultiGauss(means,covM);
}


void Application::loadTrainingSet(String dir)
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
  for(int i=0 ; i<numFiles ; ++i)
    trainingSet.push_back(loader.load(files[i]));
}

