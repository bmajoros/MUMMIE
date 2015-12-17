/**************************************************************
 baum-welch.C : For training an HMM using Expectation Maximization
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
#include "BOOM/Regex.H"
#include "BOOM/Map.H"
#include "BaumWelch.H"
#include "BaumWelchMT.H"
#include "EmissionLoader.H"
#include "SequenceSet.H"
#include "TieProfile.H"
using namespace std;

class Application {
  SequenceSet trainingSet;
  Regex filenameRegex;
  void loadMeans(const String &filename,BOOM::Vector<GSL::Vector> &);
public:
  Application();
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
  catch(int x) {
    cerr<<"exception #"<<x<<endl;
    return -1;
  }
  catch(...)
    {
      cerr << "unknown exception caught in main()" << endl;
      return -1;
    }
  return 0;
}



Application::Application()
  : filenameRegex("([^\\\\/]+)\\\\.fastb")
{
}



int Application::go(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"bB:c:C:dgIl:L:m:MN:n:rRs:t:uw:");
  if(cmd.numArgs()!=5) 
    throw "\n\
baum-welch [options] <initial.hmm> <dependency-graph.tgf> <training-dir> <#iterations> <out.hmm>\n\
   where: -b = don't \"back off\" higher-order discrete emissions\n\
          -B 0.3 = apply \"Bilmes sparsness\", keep 30% of entries\n\
          -c n = use n CPUs\n\
          -C file = center & normalize data, save transforms to file\n\
          -d = use diagonal covariance matrices\n\
          -g = use one global cov matrix\n\
          -I = use identity matrix and single coef for cov matrix\n\
          -l filename = write log-likelihood curve to file\n\
          -L threshold = stop when change in log-likelihood < threshold\n\
          -m file = initialize means from cluster file\n\
          -M = do not re-estimate means (requires -m)\n\
          -n V = add gaussian noise with variance V\n\
          -N n = use at most n training examples\n\
          -r = use one global correlation matrix (re-estimate variances)\n\
          -R = do not randomize initial HMM\n\
          -s seed = seed randomizer with given seed value\n\
          -t file = tie parameters according to profile in file\n\
          -u = update the discrete emission chains during EM\n\
          -w file = use sequence weights from file\n\
";
  String structureFile=cmd.arg(0);
  String graphFile=cmd.arg(1);
  String dir=cmd.arg(2);
  int maxIterations=cmd.arg(3).asInt();
  String outfile=cmd.arg(4);
  SparseGraph *graph=new SparseGraph(graphFile);
  int numThreads=cmd.option('c') ? cmd.optParm('c').asInt() : 1;
  ofstream *osLog=
    cmd.option('l') ? new ofstream(cmd.optParm('l').c_str()) : NULL;
  bool wantRandomize=!cmd.option('R');
  bool diagonalOnly=cmd.option('d');
  const String meansFile=cmd.option('m') ? cmd.optParm('m') : "";
  bool constantMeans=cmd.option('M');
  //if(constantMeans && !cmd.option('m')) throw "-M require -m";
  if(numThreads<0 || numThreads>20) throw "Invalid number of threads";
  const String transformFile=cmd.option('C') ? cmd.optParm('C') : "";
  double bilmesFactor=cmd.option('B') ? cmd.optParm('B').asFloat() : 1.0;
  bool useGlobalCov=cmd.option('g');
  bool useGlobalCor=cmd.option('r');
  bool useIdentityCov=cmd.option('I');
  bool wantNoise=cmd.option('n');
  double noiseVar=wantNoise ? cmd.optParm('n').asDouble() : 0;
  TieProfile *ties=cmd.option('t') ? new TieProfile(cmd.optParm('t')) : NULL;
  double LLthreshold=cmd.option('L') ? cmd.optParm('L').asFloat() : 1.0;
  bool updateDiscrete=cmd.option('u');
  bool wantBackOff=!cmd.option('b');
  String weightsFile;
  if(cmd.option('w')) weightsFile=cmd.optParm('w');

  // Misc. initialization
  cout.precision(10);
  cerr.precision(10);
  if(cmd.option('s')) SeedRandomizer(cmd.optParm('s').asUnsigned());
  else randomize();
  BOOM::catchFloatOverflow();

  // Load HMM structure
  cout<<"loading initial HMM"<<endl;
  HMM hmm(structureFile);

  // Load the training set
  cout<<"loading training data"<<endl;
  const int maxSeqs=cmd.option('N') ? cmd.optParm('N').asInt() : -1;
  trainingSet.load(dir,hmm.getSchema(),maxSeqs);
  Map<String,int> seqIDs;
  int numSeqs=trainingSet.size();
  for(int i=0 ;i<numSeqs ; ++i) {
    String filename=trainingSet[i]->getFilename();
    if(!filenameRegex.search(filename)) throw filename+" : can't parse";
    String id=filenameRegex[1];
    seqIDs[id]=i;
  }

  // Rescale the training data
  if(!transformFile.isEmpty()) {
    cout<<"rescaling data..."<<endl;
    ofstream os(transformFile.c_str());
    trainingSet.normalize();
    trainingSet.saveTransforms(os);
  }
  trainingSet.recode(hmm.getOrder()+1);

  // Add noise
  if(wantNoise) trainingSet.addNoise(0,noiseVar);

  // Load initial state
  BOOM::Vector<GSL::Vector> means;
  if(!meansFile.isEmpty()) {
    cout<<"loading initial state..."<<endl;
    loadMeans(meansFile,means);
    wantRandomize=false;
  }

  // Load sequence weights
  Array1D<float> weights(numSeqs);
  weights.setAllTo(1);
  if(weightsFile!="") {
    ifstream is(weightsFile.c_str());
    String line;
    while(!is.eof()) {
      line.getline(is);
      line.trimWhitespace();
      BOOM::Vector<BOOM::String> &fields=*line.getFields();
      if(fields.size()==2) {
	String id=fields[0];
	float weight=fields[1].asFloat();
	int seqIndex=seqIDs[id];
	weights[seqIndex]=weight;
      }
      delete &fields;
    }
  }

  // Train the HMM using the Baum-Welch algorithm
  cout << "Training..." << endl;
  BaumWelchMT bw(hmm,numThreads,maxIterations,LLthreshold,trainingSet,
		 weights,means,bilmesFactor,graph,1000000,osLog,wantRandomize,
		 diagonalOnly,useGlobalCov,useGlobalCor,useIdentityCov,
		 constantMeans,outfile,ties,updateDiscrete,wantBackOff);
  cout<<"done"<<endl;

  // Save the results
  hmm.normalizeTransitions();
  cout<<"saving hmm..."<<endl;
  hmm.save(outfile);
  cout<<"saved."<<endl;
  return 0;
}


void Application::loadMeans(const String &filename,
			    BOOM::Vector<GSL::Vector> &means)
{
  ifstream is(filename.c_str());
  if(!is.good()) throw filename+" : can't open";
  int K;
  is>>K;
  GSL::Vector v;
  for(int i=0 ; i<K ; ++i) {
    is>>v;
    means.push_back(v);
  }
}



