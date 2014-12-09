/****************************************************************
 ParallelBaumWelch.H
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/

#ifndef INCL_ParallelBaumWelch_H
#define INCL_ParallelBaumWelch_H

#include "ForwardAlgorithm.H"
#include "BackwardAlgorithm.H"
#include "HMM.H"
#include <math.h>
#include <iostream>
#include <strstream>
#include <fstream>
#include <float.h>
#include "BOOM/Vector.H"
#include "BOOM/Array1D.H"
#include "BOOM/Sequence.H"
#include "BOOM/Progress.H"
#include "BOOM/Thread.H"
using namespace std;


struct BWThread : public Thread
{
    BWThread(int states,int numAlpha,HMM &hmm,
             Vector<Sequence*> &trainingSet,int first, int last);
    virtual void f();
    void reviseExpectedEmitCounts(int sequenceLength,
                                  ForwardAlgorithm &,BackwardAlgorithm &,
                                  Sequence &);
    void reviseExpectedTransCounts(int sequenceLength,ForwardAlgorithm &,
                                   BackwardAlgorithm &,
                                   Sequence &);

    Array2D<double> A; // expected transition counts
    Array2D<double> E; // expected emission counts
    int first, last;
    double likelihoodCorrection;
    HMM &hmm;
    Vector<Sequence*> &trainingSet;
    Array1D<double> likelihoods, scalingRatios;
    int states, alphabetSize;
};


class ParallelBaumWelch
{
public:
    ParallelBaumWelch(HMM &,long maxIterations,
                      Vector<Sequence*> &trainingSet,int numThreads);
private:
    Alphabet &alphabet;
    Progress progress;
    Array1D<double> likelihoods, scalingRatios;
    Array2D<double> A; // expected transition counts
    Array2D<double> E; // expected emission counts
    HMM &hmm;
    long maxIterations;
    double logLikelihoodThreshold, sum;
    Vector<Sequence*> &trainingSet;
    int numHmmStates, numTrain, alphabetSize, numThreads;
    
    void mainAlgorithm();
};


#endif
