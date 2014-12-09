/****************************************************************
 IntervalDecoder.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "IntervalDecoder.H"
#include "BOOM/Constants.H"
#include "BOOM/SumLogProbs.H"
using namespace std;
using namespace BOOM;


IntervalDecoder::IntervalDecoder(HMM &hmm)
  : hmm(hmm), hmmGraph(hmm), numStates(hmm.countStates()),
    isFg(hmm.countStates()), FG(hmm.getForegroundStates())
{
  isFg.setAllTo(false);
  Set<int>::iterator cur=FG.begin(), end=FG.end();
  for(; cur!=end ; ++cur) isFg[*cur]=true;
}



BOOM::Vector<GffFeature> *IntervalDecoder::decode(EmissionSequence &seq,
				 const String &substrateID, 
				 const String &sourceID, 
				 const String &featureType)
{
  int L=seq.length();
  //cout<<"len="<<L<<endl;
  activeAtPos.resize(L);
  fw=new ForwardAlgorithm(hmm,seq);
  bw=new BackwardAlgorithm(hmm,seq);
  ForwardAlgorithm &F=*fw; BackwardAlgorithm &B=*bw;
  double LL=B(0,0);
  BOOM::Vector<GffFeature> *features=new BOOM::Vector<GffFeature>;
  
  // Compile set of states that have nonzero posteriors at each position
  for(int pos=0 ; pos<L ; ++pos) {
    //cout<<"pos="<<pos<<endl;
    for(int q=1 ; q<numStates ; ++q) {
      if(isFinite(F(q,pos+1)) && isFinite(B(q,pos+1))) {
	activeAtPos[pos].insert(q);
	//cout<<"  q="<<q<<endl;
      }
    }
  }

  // Iterate through all interval starting positions
  for(int begin=0 ; begin<L ; ++begin) {
    // Get set of opening states at this position
    Set<StateDoublePair> M;
    Set<int> &activeHere=activeAtPos[begin];
    for(Set<int>::iterator cur=FG.begin() ; cur!=FG.end() ; ++cur) {
      int q=*cur;
      if(!activeHere.isMember(q)) continue;
      if(begin==0) M.insert(StateDoublePair(q,F(q,begin+1),0));
      else {
	BOOM::Vector<double> V;
	const BOOM::Vector<StateDoublePair> &pred=hmmGraph.statesPreceding(q);
	for(BOOM::Vector<StateDoublePair>::const_iterator cur=pred.begin() ;
	    cur!=pred.end() ; ++cur) {
	  int p=(*cur).state;
	  if(FG.isMember(p) || !activeAtPos[begin-1].isMember(p)) continue;
	  double transP=(*cur).logP;
	  double inductiveP=F(p,begin);
	  double emitP=hmm.getEmissionProb(q,seq[begin]);
	  V.push_back(inductiveP+transP+emitP);
	}
	M.insert(StateDoublePair(q,sumLogProbs(V),0));
      }
    }
    if(M.size()==0) continue;

    // Propagate along sequence to find potential interval ends
    for(pos=begin+1 ; pos<L && M.size()>0 ; ++pos) {
      BOOM::Vector<double> V;
      for(Set<StateDoublePair>::iterator cur=M.begin();cur!=M.end();++cur) {
	StateDoublePair &sdp=(StateDoublePair&)*cur;
	if(pos==L-1) V.push_back(sdp.logP+B(sdp.state,L));
	else {
	  const BOOM::Vector<StateDoublePair> &foll=
	    hmmGraph.statesFollowing(sdp.state);
	  for(BOOM::Vector<StateDoublePair>::const_iterator cur=foll.begin() ;
	      cur!=foll.end() ; ++cur) {
	    int to=(*cur).state;
	    if(!FG.isMember(to) && activeAtPos[pos].isMember(to))
	      V.push_back(sdp.logP+(*cur).logP+
			  hmm.getEmissionProb(to,seq[pos])+B(to,pos+1)-LL);
	  }
	}
      }
      if(V.size()>0) {
	double intervalScore=exp(sumLogProbs(V));
	if(intervalScore>0) 
	  features->push_back(GffFeature("",substrateID,sourceID,featureType,
					 begin,pos,intervalScore,true,'+',0,
					 false));
      }
      Set<StateDoublePair> nextM;
      for(Set<int>::iterator cur=FG.begin() ; cur!=FG.end() ; ++cur) {
	int q=*cur;
	V.clear();
	if(!activeAtPos[pos].isMember(q)) continue;
	for(Set<StateDoublePair>::iterator cur=M.begin() ; cur!=M.end() ; 
	    ++cur) {
	  int r=(*cur).state;
	  double t=hmm.getLogTransProb(r,q);
	  if(isFinite(t))
	    V.push_back((*cur).logP+t+hmm.getEmissionProb(q,seq[pos]));
	}
	if(V.isEmpty()) continue;
	nextM.insert(StateDoublePair(q,sumLogProbs(V),0));
      }
      M=nextM;
    }
  }

  delete fw; delete bw;
  return features;
}



double IntervalDecoder::getLogP(EmissionSequence &seq,int begin,int end,
			     ForwardAlgorithm *fw,BackwardAlgorithm *bw)
{
  int L=seq.length();
  ForwardAlgorithm &F=*fw; BackwardAlgorithm &B=*bw;
  double LL=B(0,0);
  
  // Initialize leftmost matrix column
  Set<StateDoublePair> M;
  for(Set<int>::iterator cur=FG.begin() ; cur!=FG.end() ; ++cur) {
    int q=*cur;
    if(begin==0) M.insert(StateDoublePair(q,F(q,begin+1),0));
    else {
      BOOM::Vector<double> V;
      const BOOM::Vector<StateDoublePair> &pred=hmmGraph.statesPreceding(q);
      for(BOOM::Vector<StateDoublePair>::const_iterator cur=pred.begin() ;
	  cur!=pred.end() ; ++cur) {
	int p=(*cur).state;
	if(FG.isMember(p)) continue;
	double transP=(*cur).logP;
	double inductiveP=F(p,begin);
	double emitP=hmm.getEmissionProb(q,seq[begin]);
	V.push_back(inductiveP+transP+emitP);
      }
      M.insert(StateDoublePair(q,sumLogProbs(V),0));
    }
  }
  if(M.size()==0) return NEGATIVE_INFINITY;
  
  // Propagate across interval
  for(int pos=begin+1 ; pos<end ; ++pos) {
    if(M.size()==0) return NEGATIVE_INFINITY;
    Set<StateDoublePair> nextM;
    for(Set<int>::iterator cur=FG.begin() ; cur!=FG.end() ; ++cur) {
      int q=*cur;
      BOOM::Vector<double> V;
      for(Set<StateDoublePair>::iterator cur=M.begin() ; cur!=M.end() ; 
	  ++cur) {
	int r=(*cur).state;
	double t=hmm.getLogTransProb(r,q);
	if(isFinite(t))
	  V.push_back((*cur).logP+t+hmm.getEmissionProb(q,seq[pos]));
      }
      if(V.isEmpty()) continue;
      nextM.insert(StateDoublePair(q,sumLogProbs(V),0));
    }
    M=nextM;
  }

  BOOM::Vector<double> V;
  for(Set<StateDoublePair>::iterator cur=M.begin();cur!=M.end();++cur) {
    StateDoublePair &sdp=(StateDoublePair&)*cur;
    if(pos==L-1) V.push_back(sdp.logP+B(sdp.state,L));
    else {
      const BOOM::Vector<StateDoublePair> &foll=
	hmmGraph.statesFollowing(sdp.state);
      for(BOOM::Vector<StateDoublePair>::const_iterator cur=foll.begin() ;
	  cur!=foll.end() ; ++cur) {
	int to=(*cur).state;
	if(!FG.isMember(to))
	  V.push_back(sdp.logP+(*cur).logP+
		      hmm.getEmissionProb(to,seq[pos])+B(to,pos+1)-LL);
      }
    }
  }
  double logP=V.size()==0 ? NEGATIVE_INFINITY : exp(sumLogProbs(V));
  return logP;
}

