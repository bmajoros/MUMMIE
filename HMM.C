/**************************************************************
 HMM.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#include <fstream>
#include "HMM.H"
#include "BOOM/Random.H"
#include <iostream>
#include "ForwardAlgorithm.H"
#include "HMMGraph.H"
#include "BOOM/Constants.H"
#include "BOOM/RouletteWheel.H"
#include "BOOM/HigherOrderAlphabet.H"
#include "BOOM/Set.H"
using namespace std;


HMM::HMM(int numStates)
  : numStates(numStates), emissionProb(numStates),
    transitionProb(numStates,numStates), order(0)
{
  transitionProb.setAllTo(NEGATIVE_INFINITY);
}



HMM::HMM(const HMM &other)
  : emissionProb(other.emissionProb), numStates(other.numStates), 
    transitionProb(other.transitionProb), schema(other.schema),
    discreteEmitProb(other.discreteEmitProb), order(other.order),
    chains(other.chains), orders(other.orders)
{
  // copy ctor
}



HMM::HMM(const String &filename)
  : emissionProb(0), transitionProb(0,0), order(0)
{
  load(filename);
}



HMM::HMM(istream &is)
  : transitionProb(0,0), emissionProb(0), order(0)
{
  load(is);
}



HMM::~HMM()
{
  // dtor
}



const HMM &HMM::operator=(const HMM &other)
{
  emissionProb=other.emissionProb;
  numStates=other.numStates;
  transitionProb=other.transitionProb;
  schema=other.schema;
  discreteEmitProb=other.discreteEmitProb;
  order=other.order;
  chains=other.chains;
  orders=other.orders;
}



void HMM::printOn(ostream &os) const
{
  os<<"Transitions:\n";
  for(int i=0 ; i<numStates ; ++i) {
    for(int j=0 ; j<numStates ; ++j)
      os<<i<<"->"<<j<<"="<<exp(transitionProb[i][j])<<"  ";
    os<<'\n';
  }
  os<<"\nEmissions:\n";
  for(int i=0 ; i<numStates ; ++i)
    os<<i<<": "<<emissionProb[i]<<endl;
  os << endl;
}



ostream &operator<<(ostream &os,const HMM &model)
{
  model.printOn(os);
  return os;
}



double HMM::getTransitionProb(int from,int to) const
{
  return exp(transitionProb[from][to]);
}



double HMM::getLogTransProb(int from,int to) const
{
  return transitionProb[from][to];
}



void HMM::load(const String &fname)
{
  ifstream is(fname.c_str());
  if(!is.good()) throw String("Can't open ")+fname;
  load(is);
}



void HMM::load(istream &is)
{
  order=0;

  // Load number of states and allocate arrays of that size
  is>>numStates;
  discardInput(is,"states");
  emissionProb.resize(numStates);
  transitionProb.resize(numStates,numStates);

  // Load schema
  discardInput(is,"schema");
  loadSchema(is);
  const int numDiscrete=schema.getNumDiscrete();
  discreteEmitProb.resize(numStates,numDiscrete);
  chains.resize(numDiscrete);
  orders.resize(numStates,numDiscrete);

  // Read transition  && emission probabilities
  discardInput(is,"transitions");
  is >> transitionProb;
  String line;
  for(int q=1 ; q<numStates ; ++q) {
    discardInput(is,"emissions");
    emissionProb[q].load(is);
    for(int i=0 ; i<numDiscrete ; ++i) {
      int order;
      is>>order;
      orders[q][i]=order;
      if(order>this->order) this->order=order;
      discardInput(is,"order");
      discardInput(is,"alphabet");
      Alphabet alphabet;
      alphabet.load(is);
      schema.getAlphabet(i)=alphabet;
      HigherOrderAlphabet H(alphabet,order+1);
      const int numNmers=H.getNumNmers();
      Array1D<double> &row=discreteEmitProb[q][i];
      row.resize(numNmers);
      for(NmerSymbol nmer=1 ; nmer<numNmers ; ++nmer) {
	line.getline(is);
	line.trimWhitespace();
	if(is.eof() || line.isEmpty()) break;
	BOOM::Vector<BOOM::String> &fields=*line.getFields();
	Sequence seq(fields[0],alphabet);
	NmerSymbol nmer=H.lookup(seq);
	row[nmer]=fields[1].asDouble(); // in log space
	delete &fields;
      }
    }
  }

  // Load list of foreground states (if present)
  while(!is.eof()) {
    line.getline(is);
    line.trimWhitespace();
    if(line.contains("foreground")) {
      line.getline(is);
      line.trimWhitespace();
      foregroundStates.clear();
      BOOM::Vector<BOOM::String> &fg=*line.getFields(",");
      int n=fg.size();
      for(int i=0 ; i<n ; ++i) foregroundStates+=fg[i].asInt();
      delete &fg;
    }
  }	

  // Construct chains
  for(int i=0 ; i<numDiscrete ; ++i) {
    HigherOrderAlphabet H(schema.getAlphabet(i),order+1);
    chains[i]=NmerChain(H);
  }	

  // Convert to log space
  logifyTransitions();
}



Set<int> &HMM::getForegroundStates()
{
  return foregroundStates;
}



void HMM::discardInput(istream &is,const String &match)
{
  String line;
  while(!is.eof()) {
    line.getline(is);
    if(line.contains(match)) return;
  }
  throw String("Syntax error in HMM file: expecting ")+match;
}



void HMM::normalizeTransitions()
{
  unlogifyTransitions();
  double sum;
  for(int i=0 ; i<numStates ; ++i) {
    sum=0;
    for(int j=0 ; j<numStates ; ++j)
      sum+=transitionProb[i][j];
    if(sum>0)
      for(int j=0 ; j<numStates ; ++j)
	transitionProb[i][j]/=sum;
  }
  logifyTransitions();
}



void HMM::logifyTransitions()
{
  for(int i=0 ; i<numStates ; ++i) {
    Array2D<double>::RowIn2DArray<double> row=transitionProb[i];
    for(int j=0 ; j<numStates ; ++j) {
      double &cell=row[j];
      cell=log(cell);
    }
  }
}



void HMM::unlogifyTransitions()
{
  for(int i=0 ; i<numStates ; ++i) {
    Array2D<double>::RowIn2DArray<double> row=transitionProb[i];
    for(int j=0 ; j<numStates ; ++j) {
      double &cell=row[j];
      cell=exp(cell);
    }
  }
}



bool HMM::save(const String &fname) const
{
  ofstream os(fname.c_str());
  if(!os.good()) throw String("Can't create ")+fname;
  bool success=save(os);
  return success;
}



bool HMM::save(ostream &os) const
{
  unlogifyTransitions();
  os<<numStates<<" states"<<endl;
  os<<"schema:\n";
  saveSchema(os);
  os<<endl;
  os<<"transitions:\n";
  os<<transitionProb<<endl;
  int numDiscrete=schema.getNumDiscrete();
  for(int q=1 ; q<numStates ; ++q) {
    os<<"state "<<q<<" emissions:"<<endl;
    emissionProb[q].save(os); // continuous
    os<<'\n';
    for(int i=0 ; i<numDiscrete ; ++i) {
      const int order=orders[q][i];
      os<<order<<" order"<<endl;
      Alphabet &alphabet=schema.getAlphabet(i);
      os<<"alphabet:\n";alphabet.save(os);
      HigherOrderAlphabet H(alphabet,order+1);
      const int numNmers=H.getNumNmers();
      for(NmerSymbol nmer=1 ; nmer<numNmers ; ++nmer) {
	Sequence S=H.lookup(nmer);
	S.printOn(os,alphabet);os<<"\t"<<discreteEmitProb[q][i][nmer]<<endl;
      }
    }	
  }
  os<<endl;
  logifyTransitions();
  return true;
}



bool HMM::doesTransitionExist(int from,int to) const
{
  bool exists=isFinite(transitionProb[from][to]);
  return exists;
}



void HMM::addTransition(int from,int to)
{
  do
    setTransitionProb(from,to,Random0to1());
  while(!isFinite(transitionProb[from][to]));
}



int HMM::countStates() const
{
  return numStates;
}



void HMM::setTransitionProb(int from,int to,double prob)
{
  transitionProb[from][to]=log(prob);
}



double HMM::getEmissionProb(int state,const Emission &x,
			    bool omitContinuous) const
{
  double P=omitContinuous || schema.getNumContinuous()==0 ? 0.0 :
    emissionProb[state].logDensity(x.getContinuous());
  int numDiscrete=schema.getNumDiscrete();
  for(int i=0 ; i<numDiscrete ; ++i) {
    const NmerChain &chain=chains[i];
    NmerSymbol s=x.getDiscrete(i);
    const Array1D<double> &d=discreteEmitProb[state][i];
    const int dsize=d.size();
    while(int(s)>=dsize) s=chain.getSuffix(s);
    P+=d[s];
  }
  return P;
}



int HMM::numMixtureComponents() const
{
  int N=0;
  for(int q=0 ; q<numStates ; ++q) {
    int n=emissionProb[q].getNumComponents();
    if(n>N) N=n;
  }
  return N;
}



double HMM::getLogP(EmissionSequence &e) const
{
  try
  {
    HMMGraph hmmGraph(*this);
    ForwardAlgorithm f(hmmGraph,e);
    return f.getLogP();
  }
  catch(NegativeInfinityException)
  {
      return NEGATIVE_INFINITY;
  }
}



Array2D< Array1D<double> > &HMM::getDiscreteEmitDistr()
{
  return discreteEmitProb;
}



void HMM::setEmissionDistr(int state,const GaussianMixture &d)
{
  emissionProb[state]=d;
}



EmissionSequence *HMM::sample(BOOM::Vector<STATE> *statePath,int maxLen) const
{
  // Set up roulette wheels
  int nD=schema.getNumDiscrete();
  Array1D<RouletteWheel> transWheels(numStates);
  for(STATE i=0 ; i<numStates ; ++i) {
    RouletteWheel &tWheel=transWheels[i];
    for(STATE j=0 ; j<numStates ; ++j)
      tWheel.addSector(exp(transitionProb[i][j]));
    tWheel.doneAddingSectors();
  }
  Array1D<NmerChain*> chains(nD);
  Array1D<HigherOrderAlphabet> HOA(nD);
  for(int i=0 ; i<nD ; ++i) {
    HOA[i]=HigherOrderAlphabet(schema.getAlphabet(i),order+1);
    chains[i]=HOA[i].constructChain();
  }
  Array2D< Array1D<RouletteWheel> > emitWheels(numStates,nD);
  for(STATE q=1 ; q<numStates ; ++q) {
    for(int i=0 ; i<nD ; ++i) {
      NmerChain &chain=*chains[i];
      Array1D<RouletteWheel> &wheelArray=emitWheels[q][i];
      Alphabet &alphabet=schema.getAlphabet(i);
      int radix=alphabet.size();
      HigherOrderAlphabet &H=HOA[i];
      int numNmers=H.size();
      wheelArray.resize(numNmers);
      for(int nmer=0 ; nmer<numNmers ; ++nmer) {
	RouletteWheel &wheel=wheelArray[nmer];
	for(int r=0 ; r<radix ; ++r) { // r iterations thu 0th order symbols
	  NmerSymbol to=chain.nextNmer(nmer,r); // "to" is a higher order sym
	  while(int(to)>=discreteEmitProb[q][i].size())
	    to=chain.getSuffix(to);
	  wheel.addSector(exp(discreteEmitProb[q][i][to]));
	}
	wheel.doneAddingSectors();
      }
    }
  }
  // Perform sampling
  EmissionSequence *seq=new EmissionSequence(schema);
  if(statePath) statePath->clear();
  STATE prevState=0;
  Array1D<Symbol> prevEmissions(nD);
  prevEmissions.setAllTo(0); // symbol 0 is always the empty string
  int len=0, nC=schema.getNumContinuous();
  do {
    STATE nextState=transWheels[prevState].spin();
    if(nextState!=0) {
      Emission e;
      e.resize(nD,nC);
      emissionProb[nextState].sample(e.getContinuous());
      for(int i=0 ; i<nD ; ++i) {
	NmerSymbol oldNmer=prevEmissions[i];
	BaseSymbol s=emitWheels[nextState][i][oldNmer].spin();
	NmerSymbol newNmer=chains[i]->nextNmer(oldNmer,s);
	e.getDiscrete(i)=newNmer;
	prevEmissions[i]=newNmer;
      }
      seq->append(e);
      if(statePath) statePath->push_back(nextState);
    }
    prevState=nextState;
    ++len;
    if(len>maxLen) {delete seq; return NULL;}
  }
  while(prevState!=0);

  for(int i=0 ; i<nD ; ++i) delete chains[i];

  // Return sampled sequence
  return seq;
}



EmissionSequence *HMM::sampleFromState(STATE q,int length) const
{
  // Set up roulette wheels
  int nD=schema.getNumDiscrete();
  Array1D<NmerChain*> chains(nD);
  Array1D<HigherOrderAlphabet> HOA(nD);
  for(int i=0 ; i<nD ; ++i) {
    HOA[i]=HigherOrderAlphabet(schema.getAlphabet(i),order+1);
    chains[i]=HOA[i].constructChain();
  }
  Array2D< Array1D<RouletteWheel> > emitWheels(numStates,nD);
  for(STATE q=1 ; q<numStates ; ++q){
    for(int i=0 ; i<nD ; ++i) {
      NmerChain &chain=*chains[i];
      Array1D<RouletteWheel> &wheelArray=emitWheels[q][i];
      Alphabet &alphabet=schema.getAlphabet(i);
      int radix=alphabet.size();
      HigherOrderAlphabet &H=HOA[i];
      int numNmers=H.size();
      wheelArray.resize(numNmers);
      for(int nmer=0 ; nmer<numNmers ; ++nmer) {
	RouletteWheel &wheel=wheelArray[nmer];
	for(int r=0 ; r<radix ; ++r) { // r iterations thu 0th order symbols
	  NmerSymbol to=chain.nextNmer(nmer,r); // "to" is a higher order sym
	  while(int(to)>=discreteEmitProb[q][i].size())
	    to=chain.getSuffix(to);
	  wheel.addSector(exp(discreteEmitProb[q][i][to]));
	}
	wheel.doneAddingSectors();
      }
    }
  }

  // Perform sampling
  EmissionSequence *seq=new EmissionSequence(schema);
  Array1D<Symbol> prevEmissions(nD);
  prevEmissions.setAllTo(0); // symbol 0 is always the empty string
  int len=0, nC=schema.getNumContinuous();
  for(int j=0 ; j<length ; ++j) {
    Emission e;
    e.resize(nD,nC);
    emissionProb[q].sample(e.getContinuous());
    for(int i=0 ; i<nD ; ++i) {
      NmerSymbol oldNmer=prevEmissions[i];
      BaseSymbol s=emitWheels[q][i][oldNmer].spin();
      NmerSymbol newNmer=chains[i]->nextNmer(oldNmer,s);
      e.getDiscrete(i)=newNmer;
      prevEmissions[i]=newNmer;
    }
    seq->append(e);
    ++len;
  }
  for(int i=0 ; i<nD ; ++i) delete chains[i];

  // Return sampled sequence
  return seq;  
}



void HMM::saveSchema(ostream &os) const
{
  int d=schema.getNumDiscrete(), c=schema.getNumContinuous();
  os<<d<<"\t"<<c<<endl;
  for(int i=0 ; i<d ; ++i) {
    os<<schema.getDiscreteName(i)<<endl;
    schema.getAlphabet(i).save(os);os<<endl;
  }
  for(int i=0 ; i<c ; ++i) os<<schema.getContinuousName(i)<<endl;
}



void HMM::loadSchema(istream &is)
{
  int d, c;
  is>>d>>c;
  schema.setNumDiscrete(d);
  schema.setNumContinuous(c);
  String line, name;
  for(int i=0 ; i<d ; ++i) {
    is>>name;
    schema.setDiscreteName(i,name);
    line.getline(is);
    if(!schema.getAlphabet(i).load(is)) INTERNAL_ERROR;
  }
  for(int i=0 ; i<c ; ++i) {
    is>>name;
    schema.setContinuousName(i,name);
  }
}



void HMM::normalizeJoint()
{
  int numDiscrete=schema.getNumDiscrete();
  for(int i=0 ; i<numDiscrete ; ++i) {
    Alphabet &alpha=schema.getAlphabet(i);
    HigherOrderAlphabet H(alpha,order+1);
    const int numNmers=H.getNumNmers();
    const double P=1.0/numNmers;
    for(int q=1 ; q<numStates ; ++q) {
      Array1D<double> &row=discreteEmitProb[q][i];
      for(NmerSymbol nmer=1 ; nmer<numNmers ; ++nmer) row[nmer]=P;
    }
  }
}



void HMM::dropContinuousTracks()
{
  for(int i=0 ; i<numStates ; ++i)
    emissionProb[i].resize(0);
  schema.setNumContinuous(0);
}



void HMM::dropDiscreteTracks()
{
  discreteEmitProb.resize(numStates,0);
  schema.setNumDiscrete(0);
}



double HMM::getEmissionProb(int inState,const EmissionSequence &S,
			    bool omitContinuous) const
{
  int L=S.length();
  double score=0;
  for(int i=0 ; i<L ; ++i)
    score+=getEmissionProb(inState,S[i],omitContinuous);
  return score;
}



void HMM::getEquilibrium(GSL::Vector &E)
{
  GSL::TransitionMatrix T(numStates,numStates);
  for(int i=0 ; i<numStates ; ++i)
    for(int j=0 ; j<numStates ; ++j)
      T(i,j)=exp(transitionProb(i,j));
  T.getEquilibrium(E);
}



void HMM::copyEmissionDistr(const HMM &other,int otherState,int thisState)
{
  emissionProb[thisState]=other.emissionProb[otherState];
  const int numDiscrete=schema.getNumDiscrete();
  if(orders.getFirstDim()==0) {
    orders.resize(numStates,numDiscrete);
    for(int i=0 ; i<numDiscrete ; ++i) orders[0][i]=0;
  }
  for(int i=0 ; i<numDiscrete ; ++i) {
    orders[thisState][i]=other.orders[otherState][i];
    discreteEmitProb[thisState][i]=other.discreteEmitProb[otherState][i];
  }
}



void HMM::addContinuousTrack(const String &name)
{
  for(int i=0 ; i<numStates ; ++i) emissionProb[i].addVariate();
  schema.addContinuousTrack(name);
}



void HMM::addMixtureComponent()
{
  int n=emissionProb.size();
  for(int i=0 ; i<n ; ++i) emissionProb[i].addComponent();
}



void HMM::changeOrder(int o)
{
  const int numDiscrete=schema.getNumDiscrete();
  order=o;
  orders.resize(numStates,numDiscrete);
  for(int q=1 ; q<numStates ; ++q)
    for(int i=0 ; i<numDiscrete ; ++i)
      orders[q][i]=o;
}



void HMM::dropContinuousTrack(const String &name)
{
  int id=schema.lookupContinuousID(name);
  schema.dropContinuousTrack(name);
  for(int q=1 ; q<numStates ; ++q)
    emissionProb[q].dropVariate(id);
}



void HMM::dropDiscreteTrack(const String &name)
{
  int id=schema.lookupDiscreteID(name);
  schema.dropDiscreteTrack(name);
  discreteEmitProb.deleteColumn(id);
  orders.deleteColumn(id);
}



void HMM::changeStateOrder(STATE q,int newOrder)
{
  int d=schema.getNumDiscrete();
  for(int i=0 ; i<d ; ++i) {
    orders[q][i]=newOrder;
    Alphabet &alphabet=schema.getAlphabet(i);
    HigherOrderAlphabet H(alphabet,newOrder+1);
    const int numNmers=H.getNumNmers();
    Array1D<double> &discreteSlice=discreteEmitProb[q][i];
    discreteSlice.resize(numNmers);
    const double uniformP=1.0/alphabet.size();
    for(NmerSymbol nmer=0 ; nmer<numNmers ; ++nmer)
      discreteSlice[nmer]=uniformP;
  }
  normalizeDiscreteEmit();
}



void HMM::normalizeDiscreteEmit()
{
  int numDiscrete=schema.getNumDiscrete();
  for(int i=0 ; i<numDiscrete ; ++i) {
    Alphabet &alpha=schema.getAlphabet(i);
    for(int q=1 ; q<numStates ; ++q) {
      HigherOrderAlphabet H(alpha,orders[q][i]+1);
      NmerChain &chain=*H.constructChain();
      const int numNmers=H.getNumNmers();
      Array1D<double> sums(numNmers); sums.setAllTo(0.0);
      Array1D<double> &row=discreteEmitProb[q][i];
      if(row.size()!=numNmers) throw String("HMM::normalizeDiscreteEmit() ")+row.size()+" "+numNmers;
      for(NmerSymbol nmer=1 ; nmer<numNmers ; ++nmer) {
	NmerSymbol prefix=chain.getPrefix(nmer);
	sums[prefix]+=exp(row[nmer]);
      }
      for(NmerSymbol nmer=1 ; nmer<numNmers ; ++nmer) {
	NmerSymbol prefix=chain.getPrefix(nmer);
	if(isFinite(row[nmer])) row[nmer]-=log(sums[prefix]); // ###
      }
      delete &chain;
    }
  }
}





