/****************************************************************
 BaumWelch.H
 bmajoros@duke.edu

 This is OPEN SOURCE SOFTWARE governed by the ARTISTIC LICENSE.
 ****************************************************************/
#ifndef INCL_BaumWelch_H
#define INCL_BaumWelch_H
#include <math.h>
#include <iostream>
#include <strstream>
#include <fstream>
#include <float.h>
#include "BOOM/Vector.H"
#include "BOOM/Array1D.H"
#include "BOOM/Array2D.H"
#include "BOOM/SparseGraph.H"
#include "ForwardAlgorithm.H"
#include "BackwardAlgorithm.H"
#include "HMMGraph.H"
#include "EmissionSequence.H"
using namespace std;


class BaumWelch {
public:
  BaumWelch(HMM &,long maxIterations,
	    const BOOM::Vector<EmissionSequence*> &trainingSet,
	    SparseGraph *dependencyGraph=NULL,
	    int maxSampleSize=-1);
private:
  typedef pair<VertexId,VertexId> Edge;
  typedef BOOM::Vector<Edge> EdgeVector;
  const Schema &schema;
  Array1D<double> likelihoods;
  HMM &hmm;
  HMMGraph hmmGraph;
  long maxIterations;
  double logLikelihoodThreshold, sum;
  const BOOM::Vector<EmissionSequence*> &trainingSet;
  int numHmmStates, numTrain;
  int maxSampleSize;
  GSL::Vector globalMeans;
  GSL::Matrix globalCov;
  Array2D<double> A; // expected transition counts
  SparseGraph *dependencyGraph;
  EdgeVector edges;
  void mainAlgorithm();
  double computeRho(STATE,int i,int j,ForwardAlgorithm&,BackwardAlgorithm&,
		    const Emission &e);
  void getGlobalStats();
  void sparseStats();
  bool validMean(const GSL::Vector &mu);
  bool validCov(const GSL::Matrix &M);
};


#endif
