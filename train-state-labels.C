/**************************************************************
 train-state-labels.C
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
#include "BOOM/SumLogProbs.H"
#include "BOOM/GSL/Permutation.H"
#include "HMM.H"
#include "HMMGraph.H"
#include "ForwardAlgorithm.H"
#include "BackwardAlgorithm.H"
#include "EmissionLoader.H"
using namespace std;

class Application {
public:
  int go(int argc,char *argv[]);
protected:
  double evalPerm(const Permutation &,const Array2D<double> &E);
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
  CommandLine cmd(argc,argv,"p:");
  if(cmd.numArgs()!=3) 
    throw "\n\
train-state-labels [options] <input.hmm> <training-dir> <output.hmm>\n\
  where:  -p \"0 1 3 2 5 6 4\" = apply specified state mapping\n\
\n\
";
  String infile=cmd.arg(0);
  String dir=cmd.arg(1);
  String outfile=cmd.arg(2);
  bool specifiedMapping=cmd.option('p');

  // Misc. initialization
  BOOM::catchFloatOverflow();

  // Load HMM
  HMM hmm=(infile);
  Schema &schema=hmm.getSchema();
  HMMGraph hmmGraph(hmm);
  const int numStates=hmm.countStates();

  // Perform mapping
  Array1D<int> mapping(numStates);
  if(specifiedMapping) {
    String mapString=cmd.optParm('p');
    mapString.trimWhitespace();
    BOOM::Vector<String> *fields=mapString.getFields();
    if(fields->size()!=numStates) 
      throw "Mapping string contains wrong number of states";
    for(int i=0 ; i<numStates ; ++i) mapping[i]=(*fields)[i].asInt();
    delete fields;
  }
  else {
    // Process the training set
    Array2D< BOOM::Vector<double> > counts(numStates,numStates);
    BOOM::Vector<String> files;
    File::getFileList(dir,files);
    int numFiles=files.size();
    EmissionLoader loader(schema);
    double LL=0.0;
    for(int i=0 ; i<numFiles ; ++i) {
      String filename=files[i];
      int L=filename.length();
      if(L<7) continue;
      String extension=filename.substr(L-5);
      extension.toupper();
      if(extension!="FASTB") continue;
      String base=filename.substr(0,L-6);
      String pathFile=base+".path";
      EmissionSequence *S=loader.load(dir+"/"+filename);
      BOOM::Vector<int> &path=*loader.loadPathFile(dir+"/"+pathFile);
      ForwardAlgorithm F(hmmGraph,*S);
      LL+=F.getLogP();
      BackwardAlgorithm B(hmmGraph,*S);
      double logP=F.getLogP();
      int seqLen=S->length();
      for(int pos=0 ; pos<seqLen ; ++pos) {
	for(int q=1 ; q<numStates ; ++q) {
	  double posterior=F(q,pos)+B(q,pos)-logP;
	  counts[q][path[pos]].push_back(posterior);
	}
      }
      delete S;
      delete &path;
    }
    //cout<<"LL="<<LL<<endl;
    
    // Compute expected counts
    Array2D<double> expectations(numStates,numStates);
    for(int q=1 ; q<numStates ; ++q)
      for(int r=1 ; r<numStates ; ++r)
	expectations[q][r]=sumLogProbs(counts[q][r]);
    //cout<<expectations<<endl;
    
    // Find the optimal mapping
    Permutation p(numStates), bestPerm;
    double bestScore=NEGATIVE_INFINITY;
    do {
      if(p[0]>0) break;
      const double score=evalPerm(p,expectations);
      //cout<<p<<" => "<<score<<endl;
      if(score>bestScore) { bestPerm=p; bestScore=score; }
    }
    while(++p);
    for(int i=0 ; i<numStates ; ++i) mapping[i]=bestPerm[i];
  }
  cout<<"used mapping: "<<mapping<<endl;
  
  // Apply the mapping
  HMM newHMM(numStates);
  newHMM.getSchema()=hmm.getSchema();
  for(int q=1 ; q<numStates ; ++q)
    newHMM.getEmissionDistr(q)=hmm.getEmissionDistr(mapping[q]);
  for(int q=0 ; q<numStates ; ++q)
    for(int r=0 ; r<numStates ; ++r)
      newHMM.setTransitionProb(q,r,
       hmm.getTransitionProb(mapping[q],mapping[r]));

  // Save the results
  newHMM.save(outfile);

  return 0;
}



double Application::evalPerm(const Permutation &P,const Array2D<double> &E)
{
  if(P[0]!=0) return NEGATIVE_INFINITY;
  double score=0;
  const int dim=E.getFirstDim();
  for(int i=1 ; i<dim ; ++i) score+=E[i][P[i]];
  return score;
}


