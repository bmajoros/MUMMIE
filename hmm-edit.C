/**************************************************************
 hmm-edit.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#include <string>
#include <iostream>
#include <fstream>
#include "BOOM/Random.H"
#include "BOOM/Vector.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Constants.H"
#include "BOOM/Regex.H"
#include "BOOM/HigherOrderAlphabet.H"
#include "HMM.H"
using namespace std;
using namespace BOOM;

enum ParmType {
  PT_VALUES,  // array of floating-point values
  PT_ALL,     // the keyword "all"
  PT_STRING   // string representing a DNA/RNA/protein sequence
};

struct Parm {
  Parm(ParmType t) : type(t) {}
  ParmType type;
  BOOM::Vector<float> values;
  BOOM::String str;
};

class Application {
  HMM *hmm;
  Schema schema;
  BOOM::Vector<int> allStates;
  BOOM::Regex numeric, integerRange;
  Parm nextParm(CommandLine &,int &index);
  BOOM::Vector<int> getStates(Parm);
  BOOM::Vector<int> getMixtures(Parm);
  BOOM::Vector<int> getTracks(Parm);
  void trans(CommandLine &cmd,int &index);
  void mix(CommandLine &cmd,int &index);
  void mean(CommandLine &cmd,int &index);
  void cov(CommandLine &cmd,int &index);
  void var(CommandLine &cmd,int &index);
  void trk(CommandLine &cmd,int &index);
  void atrk(CommandLine &,int &index);
  void dtrk(CommandLine &cmd,int &index);
  void ddtrk(CommandLine &cmd,int &index);
  void comp(CommandLine &cmd,int &index);
  void nmer(CommandLine &cmd,int &index);
  void ord(CommandLine &cmd,int &index);
  void states(CommandLine &cmd,int &index);
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
  catch(...)
    {
      cerr << "unknown exception caught in main()" << endl;
      return -1;
    }
  return 0;
}



Application::Application() 
  : numeric("\\d"), integerRange("\\d+-\\d+")
{
}



const char *usage="\n\
hmm-edit <in-out.hmm> <operations>\n\
  where <operations> are:  \n\
    STATES n : add n states to the model\n\
    TRANS from to value : set transition probability to value\n\
    MIX state mixture_index value : set mixture weight to value\n\
    MEAN mixture_index track_index value : set mean to value\n\
    COV mixture_index track_A track_B value : set covariance to value\n\
    VAR mixture_index track_index value : set variance to value\n\
    TRK name : add continuous track\n\
    ATRK name alphabet order : add discrete track\n\
    DTRK name : delete continuous track\n\
    DDTRK name : delete discrete track\n\
    COMP : add another mixture component (applies to all states)\n\
    NMER state track_index nmer prob : set P(last nmer base | prefix)\n\
    ORD state order : set state order (for all discrete tracks)\n\
\n\
A state or index may be an individual value, comma-separated list\n\
(no spaces), or the word ALL.\n\
\n\
";

int Application::go(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  int numArgs=cmd.numArgs();
  if(numArgs<2) throw usage;
  String hmmFile=cmd.arg(0);
  cout.precision(4);
  hmm=new HMM(hmmFile);
  int numStates=hmm->countStates();
  for(int q=0 ; q<numStates ; ++q) allStates.push_back(q);
  schema=hmm->getSchema();
  for(int i=1 ; i<numArgs ; ) {
    String op=cmd.arg(i);
    ++i;
    op.toupper();
    if(op=="TRANS") trans(cmd,i);
    else if(op=="MIX") mix(cmd,i);
    else if(op=="MEAN") mean(cmd,i);
    else if(op=="COV") cov(cmd,i);
    else if(op=="VAR") var(cmd,i);
    else if(op=="TRK") trk(cmd,i);
    else if(op=="ATRK") atrk(cmd,i);
    else if(op=="DTRK") dtrk(cmd,i);
    else if(op=="DDTRK") ddtrk(cmd,i);
    else if(op=="COMP") comp(cmd,i);
    else if(op=="NMER") nmer(cmd,i);
    else if(op=="ORD") ord(cmd,i);
    else if(op=="STATES") states(cmd,i);
    else throw op+" : unrecognized command\n";
  }
  hmm->save(hmmFile);
  return 0;
}



void Application::nmer(CommandLine &cmd,int &index)
{
  // NMER state track nmer prob : set P(last nmer base | prefix)
  Parm state=nextParm(cmd,index), track=nextParm(cmd,index), 
    nmer=nextParm(cmd,index), prob=nextParm(cmd,index);
  BOOM::Vector<int> states=getStates(state);
  BOOM::Vector<int> tracks=getTracks(track);
  if(prob.type!=PT_VALUES || prob.values.size()!=1) 
    throw "single probability expected";
  if(nmer.type==PT_VALUES) throw "nmer must be non-numeric";
  if(nmer.type==PT_ALL) nmer.str="";
  Array2D< Array1D<double> > &distr=hmm->getDiscreteEmitDistr();
  for(BOOM::Vector<int>::iterator cur=states.begin(), 
	end=states.end() ; cur!=end ; ++cur) {
    int state=*cur;
    for(BOOM::Vector<int>::iterator cur=tracks.begin(), end=tracks.end() ; 
	cur!=end ; ++cur) {
      int track=*cur;
      Array1D<double> &A=distr[state][track];
      int order=hmm->getOrder(state,track);
      Alphabet &alpha=schema.getAlphabet(track);
      HigherOrderAlphabet H(alpha,order+1);
      Sequence seq(nmer.str,alpha);
      NmerSymbol s=H.lookup(seq);
      float logProb=log(prob.values[0]);
      A[s]=logProb;
      int nmerLen=nmer.str.length();
      if(nmerLen<order+1) {
	const NmerChain &chain=hmm->getChain(track);
	int numNmers=A.size();
	for(NmerSymbol S=0 ; S<numNmers ; ++S)
	  if(chain.isSuffix(S,s)) A[S]=logProb;
      }
    }
  }
}



void Application::ord(CommandLine &cmd,int &index)
{
  // ORD state track order
  Parm state=nextParm(cmd,index), ord=nextParm(cmd,index);
  int order=ord.values[0];
  BOOM::Vector<int> states=getStates(state);
  if(ord.type!=PT_VALUES || ord.values.size()!=1) 
    throw "single order expected (integer)";
  for(BOOM::Vector<int>::iterator cur=states.begin() ; cur!=states.end() ;
      ++cur) {
    int state=*cur;
    hmm->changeStateOrder(state,order);
  }
}



void Application::trans(CommandLine &cmd,int &index)
{
  // TRANS from to value : set transition probability to value
  Parm from=nextParm(cmd,index), to=nextParm(cmd,index), 
    value=nextParm(cmd,index);
  BOOM::Vector<int> fromStates=getStates(from);
  BOOM::Vector<int> toStates=getStates(to);
  if(value.type!=PT_VALUES || value.values.size()>1) 
    throw "TRANS probability must be unique";
  for(BOOM::Vector<int>::iterator cur=fromStates.begin(), 
	end=fromStates.end() ; cur!=end ; ++cur) {
    int fromState=*cur;
    for(BOOM::Vector<int>::iterator cur=toStates.begin(), end=toStates.end() ; 
	cur!=end ; ++cur) {
      int toState=*cur;
      hmm->setTransitionProb(int(from.values[0]),int(to.values[0]),
			     value.values[0]);
    }
  }
}



void Application::mix(CommandLine &cmd,int &index)
{
  // MIX state mixture_index value : set mixture weight to value
  Parm state=nextParm(cmd,index), mixIndex=nextParm(cmd,index),
    value=nextParm(cmd,index);
  BOOM::Vector<int> states=getStates(state);
  BOOM::Vector<int> mixtures=getMixtures(mixIndex);
  if(value.type!=PT_VALUES || value.values.size()>1)
    throw "MIX command must have one value";
  float Y=value.values[0];
  for(BOOM::Vector<int>::iterator cur=states.begin(), end=states.end() ; 
      cur!=end ; ++cur) {
    int state=*cur;
    if(state==0) continue;
    for(BOOM::Vector<int>::iterator cur=mixtures.begin(), end=mixtures.end() ; 
	cur!=end ; ++cur) {
      int mixIndex=*cur;
      hmm->getEmissionDistr(state).setCoef(mixIndex,Y);
    }
  }
}



void Application::mean(CommandLine &cmd,int &index)
{
  // MEAN mixture_index track_index value : set mean to value
  Parm mixIndex=nextParm(cmd,index), trackIndex=nextParm(cmd,index),
    value=nextParm(cmd,index);
  BOOM::Vector<int> mixtures=getMixtures(mixIndex);
  BOOM::Vector<int> tracks=getTracks(trackIndex);
  BOOM::Vector<int> &states=allStates;
  if(value.type!=PT_VALUES || value.values.size()>1)
    throw "MEAN command must have one value";
  float Y=value.values[0];
  for(BOOM::Vector<int>::iterator cur=states.begin(), end=states.end() ; 
      cur!=end ; ++cur) {
    int state=*cur;
    if(state==0) continue;
    for(BOOM::Vector<int>::iterator cur=tracks.begin(), end=tracks.end() ; 
	cur!=end ; ++cur) {
      int trackIndex=*cur;
      for(BOOM::Vector<int>::iterator cur=mixtures.begin(), end=mixtures.end() ; 
	  cur!=end ; ++cur) {
	int mixIndex=*cur;
	hmm->getEmissionDistr(state).getDistr(mixIndex).
	  getMeans()[trackIndex]=Y;
      }
    }
  }
}



void Application::cov(CommandLine &cmd,int &index)
{
  // COV mixture_index track_A track_B value : set covariance to value
  Parm mixIndex=nextParm(cmd,index), trackA=nextParm(cmd,index),
    trackB=nextParm(cmd,index), value=nextParm(cmd,index);
  BOOM::Vector<int> mixtures=getMixtures(mixIndex);
  BOOM::Vector<int> tracksA=getTracks(trackA), tracksB=getTracks(trackB);
  BOOM::Vector<int> &states=allStates;
  if(value.type!=PT_VALUES || value.values.size()>1)
    throw "COV command must have one value";
  float Y=value.values[0];
  for(BOOM::Vector<int>::iterator cur=states.begin(), end=states.end() ; 
      cur!=end ; ++cur) {
    int state=*cur;
    if(state==0) continue;
    for(BOOM::Vector<int>::iterator cur=tracksA.begin(), end=tracksA.end() ; 
	cur!=end ; ++cur) {
      int trackAindex=*cur;
      for(BOOM::Vector<int>::iterator cur=tracksB.begin(), end=tracksB.end() ; 
	  cur!=end ; ++cur) {
	int trackBindex=*cur;
	if(trackBindex==trackAindex) continue;
	for(BOOM::Vector<int>::iterator cur=mixtures.begin(), end=mixtures.end() ; 
	    cur!=end ; ++cur) {
	  int mixIndex=*cur;
	  hmm->getEmissionDistr(state).getDistr(mixIndex).
	    getCov()(trackAindex,trackBindex)=Y;
	}
      }
    }
  }
}



void Application::var(CommandLine &cmd,int &index)
{
  // VAR mixture_index track_index value : set variance to value
  Parm mixIndex=nextParm(cmd,index), trackIndex=nextParm(cmd,index),
    value=nextParm(cmd,index);
  BOOM::Vector<int> mixtures=getMixtures(mixIndex);
  BOOM::Vector<int> tracks=getTracks(trackIndex);
  BOOM::Vector<int> &states=allStates;
  if(value.type!=PT_VALUES || value.values.size()>1)
    throw "VAR command must have one value";
  float Y=value.values[0];
  for(BOOM::Vector<int>::iterator cur=states.begin(), end=states.end() ; 
      cur!=end ; ++cur) {
    int state=*cur;
    if(state==0) continue;
    for(BOOM::Vector<int>::iterator cur=tracks.begin(), end=tracks.end() ; 
	cur!=end ; ++cur) {
      int trackIndex=*cur;
      for(BOOM::Vector<int>::iterator cur=mixtures.begin(), end=mixtures.end() ; 
	  cur!=end ; ++cur) {
	int mixIndex=*cur;
	if(mixIndex<0 || 
	   mixIndex>=hmm->getEmissionDistr(state).getNumComponents())
	  throw "Invalid mixture index";
	hmm->getEmissionDistr(state).getDistr(mixIndex).
	  getCov()(trackIndex,trackIndex)=Y;
      }
    }
  }
}



void Application::trk(CommandLine &cmd,int &index)
{
  if(cmd.numArgs()<=index) throw "missing arguments";
  String arg=cmd.arg(index++);
  hmm->addContinuousTrack(arg);
}



void Application::atrk(CommandLine &cmd,int &index)
{
  Parm nameParm=nextParm(cmd,index), alphaParm=nextParm(cmd,index),
    orderParm=nextParm(cmd,index);
  const String name=nameParm.str;
  const String alpha=alphaParm.str;
  const int order=orderParm.str.asInt();
  hmm->addDiscreteTrack(name,Alphabet(alpha.c_str()),order);
}



void Application::dtrk(CommandLine &cmd,int &index)
{
  if(cmd.numArgs()<=index) throw "missing arguments";
  String arg=cmd.arg(index++);
  hmm->dropContinuousTrack(arg);
}



void Application::ddtrk(CommandLine &cmd,int &index)
{
  if(cmd.numArgs()<=index) throw "missing arguments";
  String arg=cmd.arg(index++);
  hmm->dropDiscreteTrack(arg);
}



void Application::comp(CommandLine &cmd,int &index)
{
  hmm->addMixtureComponent();
}



void Application::states(CommandLine &cmd,int &index)
{
  Parm parm=nextParm(cmd,index);
  const int N=parm.str.asInt();
  hmm->addStates(N);
}



Parm Application::nextParm(CommandLine &cmd,int &index)
{
  if(cmd.numArgs()<=index) throw "missing arguments";
  String arg=cmd.arg(index++);
  arg.toupper();
  arg.trimWhitespace();
  ParmType t;
  if(arg=="ALL") t=PT_ALL;
  else if(numeric.search(arg)) t=PT_VALUES;
  else t=PT_STRING;
  Parm parm(t);
  if(t==PT_VALUES) {
    BOOM::Vector<String> &fields=*arg.getFields(",");
    int numFields=fields.size();
    for(int i=0 ; i<numFields ; ++i) {
      if(integerRange.match(fields[i])) {
	int begin=integerRange[1], end=integerRange[2];
	for(int x=begin ; x<=end ; ++x)
	  parm.values.push_back(float(x));
      }
      else parm.values.push_back(fields[i].asFloat());
    }
    delete &fields;
  }
  //else if(t==PT_STRING) parm.str=arg;
  parm.str=arg;
  return parm;
}



BOOM::Vector<int> Application::getStates(Parm p)
{
  BOOM::Vector<int> v;
  if(p.type==PT_ALL) {
    int numStates=hmm->countStates();
    for(int q=0 ; q<numStates ; ++q)
      v.push_back(q);
  }
  else if(p.type==PT_VALUES) {
    int n=p.values.size();
    for(int i=0 ; i<n ; ++i)
      v.push_back(int(p.values[i]));
  }
  else throw "syntax error";
  return v;
}



BOOM::Vector<int> Application::getMixtures(Parm p)
{
  BOOM::Vector<int> v;
  if(p.type==PT_ALL) {
    int m=hmm->numMixtureComponents();
    for(int i=0 ; i<m ; ++i)
      v.push_back(i);
  }
  else if(p.type==PT_VALUES) {
    int n=p.values.size();
    for(int i=0 ; i<n ; ++i)
      v.push_back(int(p.values[i]));
  }
  else throw "syntax error";
  return v;
}



BOOM::Vector<int> Application::getTracks(Parm p)
{
  BOOM::Vector<int> v;
  if(p.type==PT_ALL) {
    int m=hmm->getSchema().getNumContinuous();
    for(int i=0 ; i<m ; ++i)
      v.push_back(i);
  }
  else if(p.type==PT_VALUES) {
    int n=p.values.size();
    for(int i=0 ; i<n ; ++i)
      v.push_back(int(p.values[i]));
  }
  else throw "syntax error";
  return v;
}


