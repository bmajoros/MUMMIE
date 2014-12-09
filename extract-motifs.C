/**************************************************************
 extract-motifs.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#include <string>
#include <iostream>
#include <fstream>
#include "BOOM/CommandLine.H"
#include "BOOM/VectorSorter.H"
#include "BOOM/Constants.H"
#include "BOOM/HigherOrderAlphabet.H"
#include "BOOM/NthOrderStringIterator.H"
#include "BOOM/Regex.H"
#include "EmissionLoader.H"
#include "SequenceSet.H"
#include "HMM.H"
using namespace std;


struct ScoredNmer {
  Sequence nmer;
  double score;
  ScoredNmer(const Sequence &nmer,double s) : nmer(nmer), score(s) {}
};

struct NmerComparator :  public Comparator<ScoredNmer> {
  bool equal(ScoredNmer &a,ScoredNmer &b)
  {return a.score==b.score;}
  bool greater(ScoredNmer &a,ScoredNmer &b)
  {return a.score>b.score;}
  bool less(ScoredNmer &a,ScoredNmer &b)
  {return a.score<b.score;}
};

class Application {
public:
  Application();
  int go(int argc,char *argv[]);
protected:
  Regex scoredRegex, unscoredRegex;
  void autoDetect(HMM &,int &fgState,int &bgState);
  void autoDetect2(HMM &,int &fgState,int &bgState);
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



Application::Application()
  : scoredRegex("(\\S+)\\s+(\\S+)"),
    unscoredRegex("(\\S+)")
{}



int Application::go(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"as:L:il:");
  if(cmd.numArgs()!=3)
    throw "\n\
extract-motifs [options] <trained.hmm> <num-motifs> <motif-length>\n\
   where: -s X = state X only\n\
          -L B = use state B as a background, compute LL ratios\n\
          -a = auto-detect foreground state, use LL ratio\n\
          -i = increment state number along length of motif\n\
          -l F = score strings from list file F\n\
\n\
";
  String hmmFile=cmd.arg(0);
  int numMotifs=cmd.arg(1).asInt();
  int motifLength=cmd.arg(2).asInt();
  bool restrictState=cmd.option('s');
  int wantState=restrictState ? cmd.optParm('s').asInt() : 0;
  int bgState=cmd.option('L') ? cmd.optParm('L').asInt() : 0;
  bool wantAuto=cmd.option('a');
  bool incrementState=cmd.option('i');
  bool wantListFile=cmd.option('l');
  String listFile;
  if(wantListFile) listFile=cmd.optParm('l');

  // Misc. initialization
  cout.precision(10);
  cerr.precision(10);
  BOOM::catchFloatOverflow();

  // Load HMM structure
  HMM hmm(hmmFile);
  const int order=hmm.getOrder();
  const int N=order+1;
  Alphabet &alpha=hmm.getSchema().getAlphabet(0);
  const int numStates=hmm.countStates();
  Schema &schema=hmm.getSchema();
  const int numTracks=schema.getNumDiscrete();
  const int numCont=schema.getNumContinuous();

  // Auto-detect
  if(wantAuto) {
    if(numStates!=3) throw "-a requires a 3-state HMM";
    if(cmd.option('L')) throw "-a and -L are mutually exclusive";
    if(cmd.option('s')) throw "-a and -s are mutually exclusive";    
    restrictState=true;
    autoDetect2(hmm,wantState,bgState);
  }

  // Some initialization
  Emission emission;
  emission.resize(numTracks,numCont);
  for(int i=0 ; i<numCont ; ++i)
    emission.getContinuous()[i]=0;

  // Load list file
  BOOM::Vector<ScoredNmer> listSeqs;
  if(wantListFile) {
    File f(listFile);
    while(!f.eof()) {
      String line=f.getline();
      if(scoredRegex.search(line)) {
	String nmerStr=scoredRegex[1];
	Sequence nmerSeq(nmerStr,alpha);
	float score=scoredRegex[2].asFloat();
	ScoredNmer scoredNmer(nmerSeq,score);
	listSeqs.push_back(scoredNmer);
      }
      else if(unscoredRegex.search(line)) {
	String nmerStr=unscoredRegex[1];
	Sequence nmerSeq(nmerStr,alpha);
	float score=scoredRegex[2].asFloat();
	ScoredNmer scoredNmer(nmerSeq,1.0);
	listSeqs.push_back(scoredNmer);
      }
    }
  }

  // Iterate through tracks
  double listScore=0;
  for(int i=0 ; i<numTracks ; ++i) {
    cout<<"track #"<<i<<endl;
    Alphabet alpha=schema.getAlphabet(i);
    HigherOrderAlphabet H(alpha,N);
    NthOrderStringIterator iter(motifLength,alpha);
    for(int q=1 ; q<numStates ; ++q) {
      if(restrictState && q!=wantState) continue;
      cout<<"state #"<<q<<endl;
      BOOM::Vector<ScoredNmer> V;
      if(wantListFile) {
	BOOM::Vector<ScoredNmer>::iterator cur=listSeqs.begin(), 
	  end=listSeqs.end();
	for(; cur!=end ; ++cur) {
	  ScoredNmer &scoredNmer=*cur;
	  Sequence &seq=scoredNmer.nmer;
	  double logP=0;
	  int qq=q;
	  motifLength=seq.getLength();
	  for(int pos=0 ; pos<motifLength ; ++pos) {
	    int begin=pos-N+1;
	    if(begin<0) begin=0;
	    NmerSymbol s=H.lookup(seq,begin,pos-begin+1);
	    emission.getDiscrete(i)=s;
	    logP+=hmm.getEmissionProb(qq,emission,true);
	    if(bgState) logP-=hmm.getEmissionProb(bgState,emission,true);
	    if(incrementState) ++qq;
	  }
	  double scoreDelta;
	  if(bgState) scoreDelta=logP*scoredNmer.score;
	  else scoreDelta=logP+log(scoredNmer.score);
	  V.push_back(ScoredNmer(seq,scoreDelta));
	  listScore+=scoreDelta;
	}
      }
      else {
	iter.reset();
	Sequence seq=iter.getNextSequence();
	for( ; !iter.done() ; seq=iter.getNextSequence()) {
	  double logP=0;
	  int qq=q;
	  for(int pos=0 ; pos<motifLength ; ++pos) {
	    int begin=pos-N+1;
	    if(begin<0) begin=0;
	    NmerSymbol s=H.lookup(seq,begin,pos-begin+1);
	    emission.getDiscrete(i)=s;
	    logP+=hmm.getEmissionProb(qq,emission,true);
	    if(bgState) logP-=hmm.getEmissionProb(bgState,emission,true);
	    if(incrementState) ++qq;
	  }
	  V.push_back(ScoredNmer(seq,logP));
	}
      }
      NmerComparator cmp;
      VectorSorter<ScoredNmer> sorter(V,cmp);
      sorter.sortDescendInPlace();
      for(int i=0 ; i<numMotifs && i<V.size() ; ++i) {
	ScoredNmer sm=V[i];
	sm.nmer.printOn(cout,alpha);
	cout<<"\t"<<sm.score<<endl;
      }
    }
  }
  if(wantListFile) cout<<"list score: "<<listScore<<endl;
  return 0;
}



void Application::autoDetect(HMM &hmm,int &fgState,int &bgState)
{
  const int numStates=hmm.countStates();
  Schema &schema=hmm.getSchema();
  const int numCont=schema.getNumContinuous();
  double bestScore=NEGATIVE_INFINITY;
  int bestState=0;
  for(int q=1 ; q<numStates ; ++q) {
    double score=0;
    GaussianMixture &mix=hmm.getEmissionDistr(q);
    int numComp=mix.getNumComponents();
    for(int j=0 ; j<numComp ; ++j) {
      double lambda=mix.getCoef(j);
      MultiGauss &gauss=mix.getDistr(j);
      const GSL::Vector &mu=gauss.getMeans();
      for(int i=0 ; i<numCont ; ++i)
	score+=lambda*mu[i];
    }
    if(score>bestScore) { bestScore=score; bestState=q; }
  }
  fgState=bestState;
  bgState=fgState==1 ? 2 : 1;
}



void Application::autoDetect2(HMM &hmm,int &fgState,int &bgState)
{
  const int numStates=hmm.countStates();
  Schema &schema=hmm.getSchema();
  const int numCont=schema.getNumContinuous();
  double bestScore=POSITIVE_INFINITY;
  int bestState=0;
  for(int q=1 ; q<numStates ; ++q) {
    double score=hmm.getTransitionProb(q,q);
    if(score<bestScore) { bestScore=score; bestState=q; }
  }
  fgState=bestState;
  bgState=fgState==1 ? 2 : 1;
}



