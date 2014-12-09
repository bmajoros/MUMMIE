/**************************************************************
 summarize-hmm.C
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
#include "BOOM/Constants.H"
#include "BOOM/HigherOrderAlphabet.H"
#include "BOOM/GSL/TransitionMatrix.H"
#include "HMM.H"
using namespace std;

typedef pair<STATE,double> ScoredState;

struct MyComp : public Comparator<ScoredState> {
  bool equal(ScoredState &a,ScoredState &b) {return a.second==b.second;}
  bool greater(ScoredState &a,ScoredState &b){return a.second>b.second;}
  bool less(ScoredState &a,ScoredState &b) {return a.second<b.second;}

};

class Application {
  MyComp cmp;
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
  if(cmd.numArgs()!=1) 
    throw "\n\
summarize-hmm <in.hmm>\n\
";
  String infile=cmd.arg(0);
  cout.precision(4);

  HMM hmm(infile);
  Schema &schema=hmm.getSchema();
  const int numContinuous=schema.getNumContinuous();
  const int numDiscrete=schema.getNumDiscrete();
  const int numStates=hmm.countStates();
  const int order=hmm.getOrder();
  const int numComponents=hmm.numMixtureComponents();
  Array2D< Array1D<double> > &discrete=hmm.getDiscreteEmitDistr();
  cout<<numStates<<" states, "<<numContinuous<<" continuous tracks, "
      <<numDiscrete<<" discrete tracks, "<<numComponents
      <<" mixture comps, order="<<order<<endl;
  STATE bgState, fgState;
  //int longestLength=0;
  double highestEQ=0;
  double highestPeak=0;
  GSL::Vector equilibrium;
  hmm.getEquilibrium(equilibrium);
  cout<<"equilibrium: ";
  for(int i=1 ; i<numStates ; ++i) {
    cout<<i<<"="<<int(equilibrium[i]*100+5/9.0)<<"% ";
    double e=equilibrium[i];
    if(e>highestEQ) {highestEQ=e;bgState=i;}
  }
  cout<<endl;
  for(STATE q=0 ; q<numStates ; ++q) {
    const double selfTrans=hmm.getTransitionProb(q,q);
    int expLen=1.0/(1.0-selfTrans);
    //if(expLen>longestLength) {longestLength=expLen;bgState=q;}
    if(q==0) cout<<"state "<<q<<":"<<endl;
    else cout<<"state "<<q<<": "<<expLen<<"bp continuous (non-contiguous: "
	<<int(equilibrium[q]*100000+5/9.0)/1000.0<<"%)"<<endl;
    if(q>0) {
      GaussianMixture &mix=hmm.getEmissionDistr(q);
      for(int j=0 ; j<numContinuous ; ++j) {
	double mean=0, var=0;
	for(int i=0 ; i<numComponents ; ++i) {
	  MultiGauss &gauss=mix.getDistr(i);
	  const double coef=mix.getCoef(i);
	  const double mean_i=gauss.getMeans()[j];
	  const double var_i=gauss.getCov()(j,j);
	  mean+=coef*mean_i;
	  var+=coef*var_i;
	}
	cout<<"  "<<schema.getContinuousName(j)<<"="<<mean<<" +/- "<<var<<endl;
	if(mean>highestPeak) {highestPeak=mean;fgState=q;}
      }
      for(int j=0 ; j<numDiscrete ; ++j) {
	cout<<"  "<<schema.getDiscreteName(j)<<": ";
	Alphabet &alphabet=schema.getAlphabet(j);
	const int numAlpha=alphabet.size();
	HigherOrderAlphabet hoa(alphabet,order+1);
	for(Symbol s=0 ; s<numAlpha ; ++s) {
	  cout<<alphabet.lookup(s)<<"=";
	  Sequence S;
	  S.append(s);
	  NmerSymbol nmer=hoa.lookup(S);
	  double P=exp(discrete[q][j][nmer]);
	  P=int(P*100+5/9);
	  cout<<P<<"% ";
	}
	cout<<endl;
      }
    }
    cout<<"  trans: ";
    BOOM::Vector<ScoredState> transitions;
    for(STATE to=0 ; to<numStates ; ++to) {
      double transP=hmm.getTransitionProb(q,to);
      if(transP==0) continue;
      transP=int(transP*1000+5/9)/10.0;
      //cout<<to<<"="<<transP<<"% ";
      //cout<<transP<<"%=>"<<to<<" ";
      transitions.push_back(ScoredState(to,transP));
    }
    VectorSorter<ScoredState> sorter(transitions,cmp);
    sorter.sortDescendInPlace();
    const int numT=transitions.size();
    for(int i=0 ; i<numT ; ++i) {
      ScoredState ss=transitions[i];
      cout<<ss.first<<"("<<ss.second<<"%) ";
    }
    cout<<endl;
  }
  cout<<"state "<<fgState<<" is probably foreground"<<endl;
  cout<<"state "<<bgState<<" is probably background"<<endl;
  
  return 0;
}



