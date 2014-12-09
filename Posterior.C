/****************************************************************
 Posterior.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "Posterior.H"
#include <iostream>
#include <math.h>
#include "HMMGraph.H"
#include "FastForward.H"
#include "BackwardAlgorithm.H"


Posterior::Posterior(HMMGraph &graph,Sequence &S,Set<int> &exonStates,
		     Set<int> &intronStates,
		     Set<int> &intergenicStates,int strandDelta)
  : dpMatrix(hmm.countStates(),S.getLength()+1),
    hmm(graph.getHMM()), 
    S(S), 
    graph(graph),
    exonStates(exonStates),
    intronStates(intronStates),
    intergenicStates(intergenicStates)
{
  // ctor

  int N=hmm.countStates();
  for(int i=0 ; i<N ; ++i)
    {
      /*
      //###
      int count=0;
      if(exonStates.isMember(i)) ++count;
      if(intronStates.isMember(i)) ++count;
      if(intergenicStates.isMember(i)) ++count;
      if(count!=1) cout<<"state "<<i<<" : "<<exonStates.isMember(i)<<" "<<intronStates.isMember(i)<<" "<<intergenicStates.isMember(i)<<endl;
      //###
      */

      int fw_state=i;
      if(fw_state>strandDelta) fw_state-=strandDelta;
      if(exonStates.isMember(fw_state)) stateToLabel[i]=0;
      else if(intronStates.isMember(fw_state)) stateToLabel[i]=1;
      else stateToLabel[i]=2;
    }
  
  compute();
}



double Posterior::operator()(int state,int pos) // uses 1-based indices!
{
  return dpMatrix[state][pos];
}



void Posterior::compute()
{
  int numStates=hmm.countStates();
  int L=S.getLength();

  // Compute forward/backward variables
  FastForward f(hmm,graph,S);
  BackwardAlgorithm b(hmm,graph,S);
  double scaledP=f.getScaledP();
  double logP=f.getLogP();

  // Compute prefix/suffix sum arrays for scaling values
  const DblArray1D &sf=f.getScalingFactors();
  const DblArray1D &sb=b.getScalingFactors();
  DblArray1D sb_sums(L+1), sf_sums(L+1);
  sf_sums[0]=0;//log(sf[0]); // ###
  sb_sums[L]=log(sb[L]);
  for(int i=1 ; i<=L ; ++i)
    sf_sums[i]=sf_sums[i-1]+log(sf[i]);
  for(int i=L-1 ; i>=0 ; --i)
    sb_sums[i]=sb_sums[i+1]+log(sb[i]);
  sb_sums[0]=sb_sums[1]; // ###

  for(int i=0 ; i<numStates ; ++i)
    {
      for(int j=1 ; j<L ; ++j)
	dpMatrix[i][j]=log(f(i,j)) + log(b(i,j)) + 
	  sf_sums[j]+sb_sums[j+1]-logP;
      dpMatrix[i][L]=log(f(i,L)) + log(b(i,L)) + 
	  sf_sums[L]-logP;
      dpMatrix[i][0]=log(0.0);
    }

  /* DEBUGGING -- ENSURE THAT THE PROBABILITIES SUM TO 1:
  for(int j=1 ; j<L ; ++j)
    {
      double sum=0;
      for(int i=0 ; i<numStates ; ++i)
	{
	  double P=exp(dpMatrix[i][j]);
	  sum+=P;
	}
      cout<<j<<" "<<sum<<endl;
    }
  */

  for(int j=1 ; j<L ; ++j)
    {
      double labelScores[3];
      labelScores[0]=labelScores[1]=labelScores[2]=0;
      
      for(int i=0 ; i<numStates ; ++i)
	{
	  double stateScore=exp(dpMatrix[i][j]);
	  labelScores[stateToLabel[i]]+=stateScore;
	}
      for(int i=0 ; i<3 ; ++i) labelScores[i]=log(labelScores[i]);
      //cout<<"exon: "<<exp(labelScores[0])<<" intron="<<exp(labelScores[1])<<" intergenic="<<exp(labelScores[2])<<endl;
      for(int i=0 ; i<numStates ; ++i)
	dpMatrix[i][j]=labelScores[stateToLabel[i]];
    }
}
