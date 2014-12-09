/****************************************************************
 BaumWelchMT.H
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_BaumWelchMT_H
#define INCL_BaumWelchMT_H
#include <math.h>
#include <iostream>
#include <strstream>
#include <fstream>
#include <float.h>
#include "BOOM/Vector.H"
#include "BOOM/Array1D.H"
#include "BOOM/Array2D.H"
#include "BOOM/Array3D.H"
#include "BOOM/SparseGraph.H"
#include "BOOM/Thread.H"
#include "BOOM/HigherOrderAlphabet.H"
#include "ForwardAlgorithm.H"
#include "BackwardAlgorithm.H"
#include "HMMGraph.H"
#include "EmissionSequence.H"
#include "SequenceSet.H"
#include "TieProfile.H"
using namespace std;

class BaumWelchMT;

/****************************************************************
                         class BaumThread1
 ****************************************************************/
struct BaumThread1 : public Thread {
  BaumThread1(const BaumWelchMT &,const BOOM::Vector<int> &seqs,
	      const int m,const int D,const int Nq,const Array1D<float> 
	      &seqWeights);
  virtual void f();
  const Array1D<float> &seqWeights;
  const BaumWelchMT &bw;
  const BOOM::Vector<int> &seqs;
  const SequenceSet &trainingSet;
  const int m, D, Nq, numDiscrete;
  const HMMGraph &hmmGraph;
  const HMM &hmm;
  DblArray1D msd;
  DblArray2D mun, ln, ld;
  Array2D<double> A;
  //Array2D< Array1D< BOOM::Vector<double> > > nmerCounts; // [state][trackNum][nmer]
  Array2D< Array1D<double> > nmerCounts; // [state][trackNum][nmer]
  double LL;
};

/****************************************************************
                         class BaumThread2
 ****************************************************************/
struct BaumThread2 : public Thread {
  BaumThread2(const BaumWelchMT &,const BOOM::Vector<int> &seqs,
	      const int m,const int D,const Array1D<float> &seqWeights);
  virtual void f();
  const Array1D<float> &seqWeights;
  const BaumWelchMT &bw;
  const HMMGraph &hmmGraph;
  const HMM &hmm;
  typedef pair<VertexId,VertexId> Edge;
  typedef BOOM::Vector<Edge> EdgeVector;
  const EdgeVector &edges;
  const BOOM::Vector<int> &seqs;
  const SequenceSet &trainingSet;
  const int m, D;
  const Array1D<GSL::Vector> &mu;
  DblArray3D sgn;
};

/****************************************************************
                        class BaumWelchMT
 ****************************************************************/
class BaumWelchMT {
public:
  BaumWelchMT(HMM &,int numThreads,long maxIterations,double LLthreshold,
	      const SequenceSet &trainingSet,const Array1D<float> &weights,
	      const BOOM::Vector<GSL::Vector> &means,double bilmesFactor=1.0,
	      SparseGraph *dependencyGraph=NULL,
	      int maxSampleSize=-1,ostream *osLog=NULL,
	      bool wantRandomize=true,bool diagonalOnly=false,
	      bool useGlobalCov=false,bool useGlobalCor=false,
	      bool useIdentityCov=false,bool constantMeans=false,
	      const String &outfile="",TieProfile* =NULL,
	      bool updateDiscrete=false,bool backOff=true);
private:
  const Array1D<float> &seqWeights;
  typedef pair<VertexId,VertexId> Edge;
  typedef BOOM::Vector<Edge> EdgeVector;
  const Schema &schema;
  TieProfile *tieProfile;
  HMM &hmm;
  String outfile;
  ostream *osLog;
  HMMGraph hmmGraph;
  long maxIterations;
  double LLthreshold, sum;
  const SequenceSet &trainingSet;
  int maxSampleSize;
  GSL::Vector globalMeans;
  GSL::Matrix globalCov; // global covariance or correlation matrix
  const BOOM::Vector<GSL::Vector> initialMeans;
  const double blimesFactor;
  Array2D<double> A; // expected transition counts
  SparseGraph *dependencyGraph;
  EdgeVector edges;
  bool wantRandomize, diagonalOnly, wantBackOff;
  bool useGlobalCov;
  bool useGlobalCor;
  bool useIdentityCov;
  bool constantMeans;
  bool wantUpdateDiscrete;
  const int numThreads;
  const int Nq;  // number of HMM states
  const int m;   // number of mixture components
  const int D;   // dimensionality of continuous emissions
  const int K;   // number of training sequences
  const int numDiscrete; // dimensionality of discrete emissions
  const double bilmesFactor; // percentage of inv cov matrix entries to keep
  DblArray1D msd;// mu & sigma denominator
  DblArray2D mun;// mu numerator
  DblArray2D ln; // lambda numerator
  DblArray2D ld; // lambda denominator
  DblArray3D sgn;// sigma numerator
  Array1D<GSL::Vector> mu; // mean vectors (one per mixture component)
  Array2D< Array1D<double> > nmerCounts; // [state][trackNum][nmer]
  Array2D<HigherOrderAlphabet> alphabets; // [trackNum]
  Array1D<NmerChain> nmerChains; // [trackNum]
  double LL;     // log likelihood
  Array1D<Thread*> threads1;
  Array1D<Thread*> threads2;
  Array1D< BOOM::Vector<int> > seqSets; // seqs used by threads
  bool shouldTieVar, shouldTieCor;
  bool shouldFixMeans, shouldFixCov;
  Array1D<bool> shouldFixWeights, shouldFixChains; // [state]
  Array2D<bool> shouldFixTransitions; // [from state][to state]
  void mainAlgorithm();
  double computeRho(STATE,int i,int j,ForwardAlgorithm&,BackwardAlgorithm&,
		    const Emission &e) const;
  void getGlobalStats();
  void getGlobalStats_fast();
  bool validMean(const GSL::Vector &mu);
  bool validCov(const GSL::Matrix &M);
  void resetCounts();
  void runThreads1();
  void runThreads2();
  void startAndWait(Array1D<Thread*> &);
  void updateCounts1();
  void updateCounts2();
  void applyBilmes(GSL::Matrix &);
  void scaledIdentity(GSL::Matrix &);
  void updateDiscreteEmit();
  void updateDiscreteEmit_joint(); // ### DEBUGGING
  void tieParms();
  void tieMeans();
  void possiblyTieMeans();
  void tieCov();
  void tieVar();
  void tieCor();
  void tieTrans(const TieTransitionEvent &);
  void tieWeights(const BOOM::Vector<STATE> &states);
  void tieChains(const BOOM::Vector<STATE> &states);
  friend class BaumThread1;
  friend class BaumThread2;
};


#endif
