/**************************************************************
 BaumWelchMT.C : the Baum-Welch Expectation Maximization algorithm
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#include <fstream>
#include "BaumWelchMT.H"
#include "BOOM/SummaryStats.H"
#include "BOOM/Correlation.H"
#include "BOOM/GSL/GaussianDistribution.H"
#include "BOOM/Array1D.H"
#include "BOOM/Array2D.H"
#include "BOOM/Array3D.H"
#include "BOOM/Constants.H"
#include "BOOM/SumLogProbs.H"
#include "BOOM/Time.H"
#include "BOOM/VectorSorter.H"


/****************************************************************
                        BaumWelchMT methods
 ****************************************************************/

BaumWelchMT::BaumWelchMT(HMM &hmm,int numThreads,long maxIterations,
			 double LLthreshold,
			 const SequenceSet &trainingSet,
			 const Array1D<float> &seqWeights,
			 const BOOM::Vector<GSL::Vector> &initialMeans,
			 double bilmesFactor,SparseGraph *dependencyGraph,
			 int maxSampleSize,ostream *osLog,
			 bool wantRandomize,bool diagonalOnly,
			 bool useGlobalCov,bool useGlobalCor,
			 bool useIdentityCov,bool constantMeans,
			 const String &outfile,TieProfile *ties,
			 bool updateDiscrete,bool backOff)
  : hmm(hmm), numThreads(numThreads), seqSets(numThreads),
    trainingSet(trainingSet),  maxIterations(maxIterations),
    hmmGraph(hmm), schema(hmm.getSchema()), maxSampleSize(maxSampleSize),
    dependencyGraph(dependencyGraph), threads1(numThreads),
    threads2(numThreads), Nq(hmm.countStates()), K(trainingSet.size()),
    D(hmm.getSchema().getNumContinuous()), m(hmm.numMixtureComponents()),
    osLog(osLog), wantRandomize(wantRandomize), diagonalOnly(diagonalOnly),
    initialMeans(initialMeans), bilmesFactor(bilmesFactor),
    useGlobalCov(useGlobalCov), useGlobalCor(useGlobalCor),
    useIdentityCov(useIdentityCov), constantMeans(constantMeans),
    outfile(outfile), numDiscrete(hmm.getSchema().getNumDiscrete()),
    tieProfile(ties), LLthreshold(LLthreshold), seqWeights(seqWeights),
    wantUpdateDiscrete(updateDiscrete), wantBackOff(backOff)
{
  // Allocate arrays
  mun.resize(D,m); 
  sgn.resize(D,D,m); 
  msd.resize(m);
  ln.resize(Nq,m); 
  ld.resize(Nq,m);
  mu.resize(m);
  globalMeans.resize(D);
  globalCov.resize(D,D);

  // Initialize alphabets
  const int numDiscrete=schema.getNumDiscrete();
  const int order=hmm.getOrder();
  const int numStates=hmm.countStates();
  nmerCounts.resize(Nq,numDiscrete);
  alphabets.resize(numStates,numDiscrete);
  nmerChains.resize(numDiscrete);
  for(int q=1 ; q<numStates ; ++q) {
    for(int i=0 ; i<numDiscrete ; ++i) {
      alphabets[q][i]=HigherOrderAlphabet(schema.getAlphabet(i),
					  hmm.getOrder(q,i)+1);
      nmerCounts[q][i].resize(alphabets[q][i].getNumNmers());
      nmerChains[i]=NmerChain(HigherOrderAlphabet(schema.getAlphabet(i),
						  hmm.getOrder()+1));
    }
  }

  // Get dependencies
  if(dependencyGraph) dependencyGraph->getEdges(edges);
  else throw "BaumWelchMT requires dependency graph";

  // Normalize emissions
  hmm.normalizeDiscreteEmit();

  // Run algorithm
  mainAlgorithm();
}



void BaumWelchMT::getGlobalStats_fast()
{
  GSL::Vector sums(D); sums.setAllTo(0);
  GSL::Matrix products(D,D); products.setAllTo(0);
  int sampleSize=0;
  for(int k=0 ; k<K ; ++k) {
    const EmissionSequence &seq=*trainingSet[k];
    const int L=seq.length();
    for(int pos=0 ; pos<L ; ++pos) {
      const Emission &s=seq[pos];
      const GSL::Vector &v=s.getContinuous();
      sums.accumulate(v);
      for(int i=0 ; i<D ; ++i)
	for(int j=i ; j<D ; ++j) {
	  //cout<<D<<" "<<i<<" "<<j<<" "<<v.getDim()<<" "<<L<<" "<<pos<<endl;
	  products(i,j)+=v[i]*v[j];
	  if(!isFinite(products(i,j))) INTERNAL_ERROR;
	}
    }
    sampleSize+=L;
  }
  for(int i=0 ; i<D ; ++i) {
    for(int j=i+1 ; j<D ; ++j)
      globalCov(i,j)=globalCov(j,i)=
	Correlation::computeR(sums[i],products(i,i),sums[j],
			      products(j,j),products(i,j),sampleSize);      
  }
  if(useGlobalCor)
    for(int i=0 ; i<D ; ++i) globalCov(i,i)=1;
  else {
    GSL::Vector sd(D);
    const double Nminus1=sampleSize-1;
    for(int i=0 ; i<D ; ++i) {
      const double sumX=sums[i];
      const double var=(products(i,i)-sumX/sampleSize*sumX)/Nminus1;
      globalCov(i,i)=var;
      sd[i]=sqrt(var);
    }
    for(int i=0 ; i<D ; ++i)
      for(int j=i+1 ; j<D ; ++j)
	globalCov(j,i)=(globalCov(i,j)*=sd[i]*sd[j]);
  }
  globalMeans=sums;
  globalMeans.scale(1.0/sampleSize);
}



void BaumWelchMT::getGlobalStats()
{
  // Compute means and variances
  globalCov.setAllTo(0.0);
  for(int d=0 ; d<D ; ++d) {
    DoubleVector V;
    for(int i=0 ; i<K ; ++i) {
      const EmissionSequence &seq=*trainingSet[i];
      const int L=seq.length();
      for(int pos=0 ; pos<L ; ++pos) {
	const Emission &s=seq[pos];
	const GSL::Vector &v=s.getContinuous();
	V.push_back(v[d]);
	if(maxSampleSize<0 || V.size()>maxSampleSize) break;
      }
    }
    SummaryStats stats(V);
    globalMeans[d]=stats.getMean();
    if(useGlobalCor) globalCov(d,d)=1.0;
    else globalCov(d,d)=stats.getVar();
  }
  
  // Compute covariances
  BOOM::Vector< pair<VertexId,VertexId> >::iterator cur=edges.begin(), 
    end=edges.end();
  for(; cur!=end ; ++cur) {
    pair<VertexId,VertexId> &p=*cur;
    VertexId d1=p.first, d2=p.second;
    if(d2==d1) continue;
    double sd1=sqrt(globalCov(d1,d1));
    double sd2=sqrt(globalCov(d2,d2));
    DoubleVector V1, V2;
    for(int i=0 ; i<K ; ++i) {
      const EmissionSequence &seq=*trainingSet[i];
      const int L=seq.length();
      for(int pos=0 ; pos<L ; ++pos) {
	const Emission &s=seq[pos];
	const GSL::Vector &v=s.getContinuous();
	V1.push_back(v[d1]);
	V2.push_back(v[d2]);
      }
      if(maxSampleSize<0 || V1.size()>maxSampleSize) break;
    }
    Correlation cor(V1,V2);
    if(useGlobalCor)
      globalCov(d2,d1)=globalCov(d1,d2)=cor.getR();
    else
      globalCov(d2,d1)=globalCov(d1,d2)=cor.getR()*sd1*sd2;
  }
}



double BaumWelchMT::computeRho(STATE q,int i,int j,ForwardAlgorithm &F,
			     BackwardAlgorithm &B,const Emission &e) const
{
  ++i;
  const GaussianMixture &mix=hmm.getEmissionDistr(q);
  const GSL::Vector &Si=e.getContinuous();
  BOOM::Vector<double> V;
  for(STATE s=1 ; s<Nq ; ++s) V.push_back(safeAdd(F(s,i),B(s,i)));
  double den=sumLogProbs(V);
  V.clear();
  for(int k=0 ; k<m ; ++k) 
    V.push_back(safeAdd(mix.getLogCoef(k),mix.getDistr(k).logDensity(Si)));
  den+=sumLogProbs(V);
  double num=
    F(q,i)+B(q,i)+mix.getLogCoef(j)+mix.getDistr(j).logDensity(Si);
  double rho=exp(num-den);
  if(!isFinite(rho) || rho<0 || rho>1.001) {
    cout<<"rho="<<rho<<endl;
    INTERNAL_ERROR;
  }
  return rho;
}



void BaumWelchMT::resetCounts()
{
  msd.setAllTo(0.0);
  mun.setAllTo(0.0);
  ln.setAllTo(0.0);
  ld.setAllTo(0.0);
  sgn.setAllTo(0.0);
  A.setAllTo(0);
  for(int q=1 ; q<Nq ; ++q)
    for(int i=0 ; i<numDiscrete ; ++i) 
      nmerCounts[q][i].setAllTo(NEGATIVE_INFINITY);
  shouldTieVar=shouldTieCor=false;
  shouldFixMeans=shouldFixCov=false;
  shouldFixWeights.resize(Nq);
  shouldFixWeights.setAllTo(false);
  shouldFixChains.resize(Nq);
  shouldFixChains.setAllTo(false);
  shouldFixTransitions.resize(Nq,Nq);
  shouldFixTransitions.setAllTo(false);
}



void BaumWelchMT::runThreads1()
{
  for(int t=0 ; t<numThreads ; ++t)
    threads1[t]=new BaumThread1(*this,seqSets[t],m,D,Nq,seqWeights);
  startAndWait(threads1);
}



void BaumWelchMT::runThreads2()
{
  for(int t=0 ; t<numThreads ; ++t)
    threads2[t]=new BaumThread2(*this,seqSets[t],m,D,seqWeights);
  startAndWait(threads2);
}



void BaumWelchMT::startAndWait(Array1D<Thread*> &threads)
{
  for(int t=0 ; t<numThreads ; ++t) threads[t]->start();
  for(int t=0 ; t<numThreads ; ++t) threads[t]->join();
}



void BaumWelchMT::updateCounts1()
{
  for(int i=0 ; i<numThreads ; ++i) {
    BaumThread1 &t=static_cast<BaumThread1&>(*threads1[i]);
    LL+=t.LL;
    for(int j=0 ; j<m ; ++j) {
      msd[j]+=t.msd[j];
      for(int d=0 ; d<D ; ++d)
	mun[d][j]+=t.mun[d][j];
      for(int q=1 ; q<Nq ; ++q) {
	ln[q][j]+=t.ln[q][j];
	ld[q][j]+=t.ld[q][j];
      }
    }
    for(int q=0 ; q<Nq ; ++q)
      for(int r=0 ; r<Nq ; ++r)
	A[q][r]+=t.A[q][r];
  }
  //nmerCounts.setAllTo(NEGATIVE_INFINITY);
  for(int q=1 ; q<Nq ; ++q) {
    for(int i=0 ; i<numDiscrete ; ++i) {
      Array1D<double> &row=nmerCounts[q][i];
      const int numNmers=row.size();
      for(int j=0 ; j<numNmers ; ++j) {
	/*
	BOOM::Vector<double> V;
	for(int k=0 ; k<numThreads ; ++k) {
	  BaumThread1 &t=static_cast<BaumThread1&>(*threads1[k]);
	  V.append(t.nmerCounts[q][i][j]);
	}
	row[j]=sumLogProbs(V);
	*/
	for(int k=0 ; k<numThreads ; ++k) {
	  BaumThread1 &t=static_cast<BaumThread1&>(*threads1[k]);
	  row[j]=sumLogProbs(row[j],t.nmerCounts[q][i][j]);
	}
      }
    }
  }
  for(int i=0 ; i<numThreads ; ++i) delete threads1[i];
  cout.precision(5);
  cout<<"RHO(";
  for(int i=0 ; i<m ; ++i) {
    cout<<int(msd[i]+5.0/9);
    if(i<m-1) cout<<" ";
  }
  cout<<")\t";
}



void BaumWelchMT::updateCounts2()
{
  for(int i=0 ; i<numThreads ; ++i) {
    BaumThread2 &t=static_cast<BaumThread2&>(*threads2[i]);
    BOOM::Vector< pair<VertexId,VertexId> >::iterator cur=edges.begin(), 
      end=edges.end();
    for(; cur!=end ; ++cur) {
      pair<VertexId,VertexId> &p=*cur;
      const VertexId a=p.first, b=p.second;
      DblArray3D::IndexedTwice<double> s=sgn[a][b], z=sgn[b][a], 
	ts=t.sgn[a][b], tz=t.sgn[b][a];
      for(int j=0 ; j<m ; ++j) { s[j]+=ts[j]; /*if(a!=b) z[j]+=tz[j]; */}
    }
    delete &t;
  }
}



void BaumWelchMT::mainAlgorithm()
{
  // Observe global stats
  cout<<"initializing global stats"<<endl;
  getGlobalStats_fast();
  Array1D<GaussianDistribution> gauss(D);//###should use MultiGauss
  for(int d=0 ; d<D ; ++d) 
    gauss[d]=GaussianDistribution(globalMeans[d],globalCov(d,d));

  // Initialize the means and covariance matrices
  cout<<"initializing means and covariance matrices..."<<endl;
  const int numTrain=trainingSet.size();
  for(int j=0 ; j<m ; ++j) {
    GSL::Vector &mu_j=mu[j];
    mu_j.resize(D);
    if(wantRandomize) {
      for(int d=0 ; d<D ; ++d) 	mu_j[d]=gauss[d].random();
      /*
      const EmissionSequence &Si=*trainingSet[RandomNumber(numTrain)];
      const int L=Si.length();
      mu_j=Si[RandomNumber(L)].getContinuous();
      */
    }
    else if(initialMeans.size()>j) mu_j=initialMeans[j];
    else mu_j=hmm.getEmissionDistr(1).getDistr(j).getMeans();
  }
  if(wantRandomize || initialMeans.size()>0)
    for(int j=0 ; j<m ; ++j) {
      MultiGauss gauss(mu[j],globalCov,diagonalOnly);
      for(STATE q=1 ; q<Nq ; ++q) {
	GaussianMixture &mix=hmm.getEmissionDistr(q);
	mix.setDistr(j,gauss);
       }
    }

  // Initialize lambda
  cout<<"initializing mixture weights..."<<endl;
  if(wantRandomize || initialMeans.size()>0)
    for(STATE q=1 ; q<Nq ; ++q) {
      GaussianMixture &mix=hmm.getEmissionDistr(q);
      Array1D<double> lambda(m);
      double sum=0;
      for(int j=0 ; j<m ; ++j) sum+=lambda[j]=Random0to1();
      for(int j=0 ; j<m ; ++j) mix.setCoef(j,lambda[j]/sum);
    }
  
  // Initialize transition counts
  A.resize(Nq,Nq);

  // Initialize sequence sets for the threads
  cout<<"initializing thread sets..."<<endl;
  float numSeqsPerThread=K/float(numThreads), first=0;
  for(int t=0 ; t<numThreads ; ++t) {
    int next=int(first+numSeqsPerThread);
    if(next>K) next=K;
    BOOM::Vector<int> &seqSet=seqSets[t];
    for(int k=int(first) ; k<next ; ++k)
      seqSet.push_back(k);
    first=next;
  }

  //cout<<"saving initial hmm..."<<endl;
  if(!outfile.isEmpty()) {
    hmm.normalizeTransitions();
    //hmm.save(outfile);
  }

  //====================================================== EM EM EM EM
  //====================================================== EM EM EM EM
  //====================================================== EM EM EM EM
  // Perform EM
  cout<<"performing EM"<<endl;
  double prevLL;
  int violations=0;
  for(int iter=0 ; iter<maxIterations ; ++iter) {
    LL=0.0;
    cout<<getDateAndTime();
    cout<<"ITERATION #"<<iter<<"\t";cout.flush();
    resetCounts();
    hmmGraph.updateProbs();

    // Run the first set of threads
    runThreads1();
    updateCounts1();

    if(tieProfile) possiblyTieMeans();

    // Update means (but don't change HMM yet)
    if(!constantMeans && !shouldFixMeans) {
      for(int j=0 ; j<m ; ++j)
	for(int d=0 ; d<D ; ++d) {
	  const double den=msd[j];
	  double mu_j=(den==0 || !isFinite(den)) ? 0 : mun[d][j]/den;
	  mu[j][d]=mu_j;
	}
    }

    // Run the second set of threads
    if(!useGlobalCov) {
      runThreads2();
      updateCounts2();
    }

    // Perform parameter tying, if requested
    if(tieProfile) tieParms();

    // Update transition probabilities
    for(int k=0 ; k<Nq ; ++k) {
      double sum=0;
      const BOOM::Vector<StateDoublePair> &next=hmmGraph.statesFollowing(k);
      BOOM::Vector<StateDoublePair>::const_iterator cur=next.begin(), 
	end=next.end();
      for(; cur!=end ; ++cur) {
	const StateDoublePair &sdp=*cur;
	STATE l=sdp.state;
	sum+=A[k][l];
      }
      cur=next.begin();
      if(sum==0) {
	for(; cur!=end ; ++cur)
	  if(!shouldFixTransitions[k][(*cur).state])
	    hmm.setTransitionProb(k,(*cur).state,0.0);
      }
      else for(; cur!=end ; ++cur) {
	const int l=(*cur).state;
	if(shouldFixTransitions[k][l]) continue;
	double p=A[k][l]/sum;
	//if(p<1e-10) p=0.0; // ###
	hmm.setTransitionProb(k,l,p);
      }
    }
    
    // Update mixture coefficients
    for(int j=0 ; j<m ; ++j) {
      const bool msd_bad=msd[j]==0 || !isFinite(msd[j]);
      for(STATE q=1 ; q<Nq ; ++q) {
	if(shouldFixWeights[q]) continue;
	const double den=ld[q][j];
	double lambda=
	  (den==0 || !isFinite(den) || msd_bad) ? 0.0 : ln[q][j]/den;
	if(!isFinite(lambda)) INTERNAL_ERROR;
	hmm.getEmissionDistr(q).setCoef(j,lambda);
      }
    }
    
    // Update covariance matrix
    //   if(!shouldFixMeans && !shouldFixCov) {
    for(int j=0 ; j<m ; ++j) {
      GSL::Matrix Cj(D,D);
      if(useGlobalCor) {
	GSL::Vector sd(D);
	const double den=msd[j];
	for(int d=0 ; d<D ; ++d)
	  sd[d]=sqrt(sgn[d][d][j]/den);
	for(int c=0 ; c<D ; ++c)
	  for(int d=0 ; d<D ; ++d)
	    Cj(c,d)=globalCov(c,d)*sd[c]*sd[d];
      }
      else if(useGlobalCov) Cj=globalCov;
      else {
	const double den=msd[j];
	Cj.setAllTo(0.0);
	EdgeVector::iterator cur=edges.begin(), end=edges.end();
	for(; cur!=end ; ++cur) {
	  const Edge &edge=*cur;
	  VertexId d=edge.first, c=edge.second;
	  Cj(d,c)=(den==0 || !isFinite(den)) ? 1 : sgn[d][c][j]/den;
	  if(c!=d) Cj(c,d)=Cj(d,c);
	}
	if(bilmesFactor!=1.0) applyBilmes(Cj);
	else if(useIdentityCov) scaledIdentity(Cj);
      }	
      GSL::Vector means=mu[j];
      GSL::Matrix cov=Cj;
      if(shouldFixMeans) 
	means=hmm.getEmissionDistr(1).getDistr(j).getMeans();
      if(shouldFixCov)
	cov=hmm.getEmissionDistr(1).getDistr(j).getCov();
      MultiGauss gauss(means,cov,diagonalOnly);
      for(STATE q=1 ; q<Nq ; ++q)
	hmm.getEmissionDistr(q).setDistr(j,gauss);
    }
    //	}
    if(shouldTieVar) tieVar();
    if(shouldTieCor) tieCor();
    if(wantUpdateDiscrete) updateDiscreteEmit();

    // Prepare for next iteration
    cout.precision(10);
    if(osLog) osLog->precision(10);
    if(osLog) (*osLog)<<iter<<"\t"<<LL<<endl;
    const double deltaLL=LL-prevLL;
    cout<<"LL="<<LL;
    if(iter>0) cout<<"  deltaLL="<<deltaLL;
    if(iter>0 && LL<prevLL) {
      cout<<"  <<<<< LL DECREASE";
      ++violations;
    }
    cout<<endl;
    if((iter>0 && deltaLL<=LLthreshold) || iter+1>=maxIterations) {
      if(!wantUpdateDiscrete) updateDiscreteEmit();
      break;
    }
    prevLL=LL;
    if(!outfile.isEmpty()) {
      hmm.normalizeTransitions();
      String snapshot=outfile;
      if(snapshot.length()>0 && !snapshot.contains("/") && 
	 File::exists("snapshot"))
	snapshot=String("snapshot/")+outfile;
      hmm.save(snapshot+"-snapshot"+iter);
    }
  }
  cout<<"total LL violations: "<<violations<<endl;
  cout<<getDateAndTime()<<endl;
}


void BaumWelchMT::updateDiscreteEmit()
{
  cout<<"updating DNA track"<<endl;
  Array2D< Array1D<double> > &distr=hmm.getDiscreteEmitDistr();
  for(int q=1 ; q<Nq ; ++q) {
    if(shouldFixChains[q]) continue;
    for(int i=0 ; i<numDiscrete ; ++i) {
      Array1D<double> &counts=nmerCounts[q][i], &prob=distr[q][i];
      HigherOrderAlphabet &H=alphabets[q][i];
      const Alphabet &A=H.getBaseAlphabet();
      const int numAlpha=A.size();
      const int numNmers=H.getNumNmers();
      NmerChain &chain=nmerChains[i];
      for(NmerSymbol nmer=1 ; nmer<numNmers ; ++nmer) {
	NmerSymbol prefix=chain.getPrefix(nmer);
	//prob[nmer]=counts[nmer]-counts[prefix];
	BOOM::Vector<double> V;
	for(BaseSymbol s=0 ; s<numAlpha ; ++s) {
	  NmerSymbol otherNmer=chain.nextNmer(prefix,s);
	  V.push_back(counts[otherNmer]);
	}
	const double denom=sumLogProbs(V);
	if(isFinite(denom)) prob[nmer]=counts[nmer]-denom;
	else if(wantBackOff) 
	  prob[nmer]=prob[chain.getSuffix(nmer)]; // "back off" ###
	else prob[nmer]=NEGATIVE_INFINITY;
      }
    }
  }
}



void BaumWelchMT::updateDiscreteEmit_joint()
{
  cout<<"updating DNA track (joint)"<<endl;
  Array2D< Array1D<double> > &distr=hmm.getDiscreteEmitDistr();
  for(int q=1 ; q<Nq ; ++q) {
    if(shouldFixChains[q]) continue;
    for(int i=0 ; i<numDiscrete ; ++i) {
      Array1D<double> &counts=nmerCounts[q][i], &prob=distr[q][i];
      HigherOrderAlphabet &H=alphabets[q][i];
      const int numNmers=H.getNumNmers();
      BOOM::Vector<double> V;
      for(NmerSymbol nmer=1 ; nmer<numNmers ; ++nmer) 
	V.push_back(counts[nmer]);
      double divisor=sumLogProbs(V);
      for(NmerSymbol nmer=1 ; nmer<numNmers ; ++nmer)
	prob[nmer]=counts[nmer]-log(divisor);
    }
  }
}



bool BaumWelchMT::validMean(const GSL::Vector &mu) 
{
  for(int i=0 ; i<D ; ++i) if(!isFinite(mu[i])) return false;
  return true;
}



bool BaumWelchMT::validCov(const GSL::Matrix &M)
{
  for(int i=0 ; i<D ; ++i)
    for(int j=0 ; j<D ; ++j)
      if(!isFinite(M(i,j))) return false;
  return true;
}



void BaumWelchMT::applyBilmes(GSL::Matrix &C)
{
  GSL::Matrix Cinv;
  C.invert(Cinv);
  DoubleVector V;
  const int n=C.getNumRows();
  for(int i=0 ; i<n ; ++i)
    for(int j=0 ; j<n ; ++j)
      V.push_back(C(i,j));
  DirectComparator<double> cmp;
  VectorSorter<double> sorter(V,cmp);
  sorter.sortDescendInPlace();
  const int nV=V.size();
  int index=blimesFactor*nV;
  if(index>=nV) index=nV-1;
  const double cutoff=V[index];
  for(int i=0 ; i<n ; ++i)
    for(int j=0 ; j<n ; ++j) {
      double &e=Cinv(i,j);
      if(e<cutoff) e=0.0;
    }
  Cinv.invert(C);
}



void BaumWelchMT::scaledIdentity(GSL::Matrix &M)
{
  // Make sure it's diagonal
  const int n=M.getNumRows();
  for(int i=0 ; i<n ; ++i)
    for(int j=0 ; j<n ; ++j)
      if(i!=j) M(i,j)=0.0;

  // Compute mean diag value
  double sum=0;
  for(int i=0 ; i<n ; ++i) sum+=M(i,i)/n;

  // Set all diagonal entries to mean
  for(int i=0 ; i<n ; ++i) M(i,i)=sum;
}



void BaumWelchMT::tieParms()
{
  TieProfile &tieProfile=*this->tieProfile;
  const int numEvents=tieProfile.numEvents();
  for(int i=0 ; i<numEvents ; ++i) {
    const TieEvent &E=tieProfile[i];
    if(E.getOp()==TIE_PARMS) {
      switch(E.getEntity()) {
      case TIE_MEANS: tieMeans(); break;
      case TIE_COV_MATRIX: tieCov(); break;
      case TIE_VARIANCES: shouldTieVar=true; break;
      case TIE_COR_MATRIX: shouldTieCor=true; break;
      case TIE_TRANSITIONS: 
	tieTrans(static_cast<const TieTransitionEvent&>(E)); break;
      case TIE_MIXTURE_WEIGHTS: tieWeights(E.getStates()); break;
      case TIE_MARKOV_CHAINS: tieChains(E.getStates()); break;
      default: INTERNAL_ERROR;
      }
    }
    else { // op==FIX
      switch(E.getEntity()) {
      case TIE_MEANS: shouldFixMeans=true; break;
      case TIE_COV_MATRIX: shouldFixCov=true; break;
      case TIE_TRANSITIONS: {
	const TieTransitionEvent &TE=static_cast<const TieTransitionEvent&>(E);
	const int n=TE.getNumTransitions();
	STATE from, to;
	for(int j=0 ; j<n ; ++j) {
	  TE.getTransition(j,from,to);
	  if(from<0 || to<0)  // means "all transitions"
	    shouldFixTransitions.setAllTo(true);
	  else shouldFixTransitions[from][to]=true;
	}} break;
      case TIE_MIXTURE_WEIGHTS: {
	const BOOM::Vector<STATE> &Q=E.getStates();
	const int n=Q.size();
	if(n==0) shouldFixWeights.setAllTo(true);
	else for(int j=0 ; j<n ; ++j) shouldFixWeights[Q[j]]=true;
      } break;
      case TIE_MARKOV_CHAINS: {
	const BOOM::Vector<STATE> &Q=E.getStates();
	const int n=Q.size();
	if(n==0) shouldFixChains.setAllTo(true);
	else for(int j=0 ; j<n ; ++j) shouldFixChains[Q[j]]=true;
      } break;
      case TIE_VARIANCES: // fall through
      case TIE_COR_MATRIX:// fall through
      default: INTERNAL_ERROR;
      }
    }
  }
}



void BaumWelchMT::possiblyTieMeans()
{
  TieProfile &tieProfile=*this->tieProfile;
  const int numEvents=tieProfile.numEvents();
  for(int i=0 ; i<numEvents ; ++i) {
    const TieEvent &E=tieProfile[i];
    if(E.getOp()!=TIE_PARMS) continue;
    if(E.getEntity()==TIE_MEANS) {
      tieMeans();
      return;
    }
  }
}



void BaumWelchMT::tieMeans()
{
  if(constantMeans) return;
  double den=0;
  for(int j=0 ; j<m ; ++j) den+=msd[j];
  den/=m;
  for(int j=0 ; j<m ; ++j) msd[j]=den;
  for(int d=0 ; d<D ; ++d) {
    double numer=0;
    for(int j=0 ; j<m ; ++j)  numer+=mun[d][j];
    numer/=m;
    for(int j=0 ; j<m ; ++j) mun[d][j]=numer;
  }
}



void BaumWelchMT::tieCov()
{
  double denSum=0;
  for(int j=0 ; j<m ; ++j) denSum+=msd[j];
  denSum/=m;
  for(int j=0 ; j<m ; ++j) msd[j]=denSum;
  for(int d=0 ; d<D ; ++d) {
    for(int e=0 ; e<D ; ++e) {
      double numerSum=0;
      for(int j=0 ; j<m ; ++j) numerSum+=sgn[d][e][j];
      numerSum/=m;
      for(int j=0 ; j<m ; ++j) sgn[d][e][j]=numerSum;
    }
  }
}



void BaumWelchMT::tieVar()
{
  GaussianMixture &mix=hmm.getEmissionDistr(1);
  const int m=mix.getNumComponents();
  GSL::Vector sdAve(D); sdAve.setAllTo(0);
  for(int j=0 ; j<m ; ++j) {
    MultiGauss &gauss=mix.getDistr(j);
    GSL::Matrix &Cov=gauss.getCov(), Cor;
    GSL::Vector sd;
    Cov.covToCor(Cor,sd);
    sdAve.accumulate(sd);
  }
  sdAve.scale(1.0/m);
  for(int j=0 ; j<m ; ++j) {
    MultiGauss &gauss=mix.getDistr(j);
    GSL::Matrix &Cov=gauss.getCov(), Cor;
    GSL::Vector sd;
    Cov.covToCor(Cor,sd);
    Cov.rebuildCovMatrix(Cor,sdAve);
  }
  for(STATE q=1 ; q<Nq ; ++q) 
    hmm.getEmissionDistr(q)=mix;
}



void BaumWelchMT::tieCor()
{
  GaussianMixture &mix=hmm.getEmissionDistr(1);
  const int m=mix.getNumComponents();
  GSL::Matrix corAve(D,D); corAve.setAllTo(0);
  for(int j=0 ; j<m ; ++j) {
    MultiGauss &gauss=mix.getDistr(j);
    GSL::Matrix &Cov=gauss.getCov(), Cor;
    GSL::Vector sd;
    Cov.covToCor(Cor,sd);
    corAve+=(Cor);
  }
  corAve.scale(1.0/m);
  for(int j=0 ; j<m ; ++j) {
    MultiGauss &gauss=mix.getDistr(j);
    GSL::Matrix &Cov=gauss.getCov(), Cor;
    GSL::Vector sd;
    Cov.covToCor(Cor,sd);
    Cov.rebuildCovMatrix(corAve,sd);
  }
  for(STATE q=1 ; q<Nq ; ++q) 
    hmm.getEmissionDistr(q)=mix;
}



void BaumWelchMT::tieTrans(const TieTransitionEvent &e)
{
  const int n=e.getNumTransitions();
  double sum=0;
  STATE from, to;
  for(int i=0 ; i<n ; ++i) {
    e.getTransition(i,from,to);
    sum+=A[from][to];
  }
  const double ave=sum/n;
  for(int i=0 ; i<n ; ++i) {
    e.getTransition(i,from,to);
    A[from][to]=ave;
  }
}



void BaumWelchMT::tieWeights(const BOOM::Vector<STATE> &states)
{
  const int n=states.size();
  for(int j=0 ; j<m ; ++j) {
    double numer=0, den=0;
    for(int i=0 ; i<n ; ++i) {
      const STATE q=states[i];
      numer+=ln[q][j];
      den+=ld[q][j];
    }
    numer/=n;
    den/=n;
    for(int i=0 ; i<n ; ++i) {
      const STATE q=states[i];
      ln[q][j]=numer;
      ld[q][j]=den;
    }
  }
}



void BaumWelchMT::tieChains(const BOOM::Vector<STATE> &states)
{
  const int numStates=states.size();
  for(int d=0 ; d<D ; ++d) {
    const int numNmers=alphabets[states[0]][d].getNumNmers();
    Array1D<double> sums(numNmers); sums.setAllTo(0);
    for(int i=0 ; i<numStates ; ++i)
      for(NmerSymbol s=0 ; s<numNmers ; ++s) 
	sums[s]+=nmerCounts[states[i]][d][s];
    for(int i=0 ; i<numStates ; ++i)
      for(NmerSymbol s=0 ; s<numNmers ; ++s) 
	nmerCounts[states[i]][d][s]=sums[s]/numStates;
  }	
}



/****************************************************************
                          BaumThread1 methods
 ****************************************************************/
BaumThread1::BaumThread1(const BaumWelchMT &parent,
			 const BOOM::Vector<int> &seqs,
			 const int m,const int D,const int Nq,
			 const Array1D<float> &seqWeights)
  : bw(parent), seqs(seqs), msd(m), mun(D,m), ln(Nq,m), ld(Nq,m),
    m(m), D(D), Nq(Nq), A(Nq,Nq), trainingSet(parent.trainingSet),
    hmmGraph(parent.hmmGraph), hmm(parent.hmm),
    numDiscrete(parent.schema.getNumDiscrete()),
    seqWeights(seqWeights)
{
  msd.setAllTo(0.0);
  mun.setAllTo(0.0);
  ln.setAllTo(0.0);
  ld.setAllTo(0.0);
  A.setAllTo(0.0);
  LL=0.0;
  nmerCounts.resize(Nq,numDiscrete);
  for(int i=0 ; i<numDiscrete ; ++i)
    for(int q=1 ; q<Nq ; ++q) {
      //Array1D< BOOM::Vector<double> > &row=nmerCounts[q][i];
      Array1D<double> &row=nmerCounts[q][i];
      row.resize(parent.alphabets[q][i].getNumNmers());
      row.setAllTo(NEGATIVE_INFINITY);
    }
}



void BaumThread1::f()
{
  int numSeqs=seqs.size();
  BOOM::Vector<int>::const_iterator cur=seqs.begin(), end=seqs.end();
  for(; cur!=end ; ++cur) {
    const int k=*cur;
    const EmissionSequence &S=*trainingSet[k];
    const float seqWeight=seqWeights[k];
    //S.save(cout);
    const int L=S.length();
    ForwardAlgorithm F(hmmGraph,S);
    BackwardAlgorithm B(hmmGraph,S);
    for(int j=0 ; j<m ; ++j) {
      for(STATE q=1 ; q<Nq ; ++q) {
	for(int i=0 ; i<L ; ++i) {
	  const Emission &Si=S[i];
	  const double rho=bw.computeRho(q,i,j,F,B,Si);
	  ln[q][j]+=rho*seqWeight;
	  msd[j]+=rho*seqWeight;
	  for(int k=0 ; k<m ; ++k) 
	    ld[q][j]+=bw.computeRho(q,i,k,F,B,Si)*seqWeight;
	  for(int d=0 ; d<D ; ++d) mun[d][j]+=rho*Si[d]*seqWeight;
	}
      }
    }//foreach component j...

    // Update transition counts
    double seqP=F.getLogP();
    LL+=seqP;
    for(int k=0 ; k<Nq ; ++k) {
      const BOOM::Vector<StateDoublePair> &next=hmmGraph.statesFollowing(k);
      BOOM::Vector<StateDoublePair>::const_iterator cur=next.begin(), 
	end=next.end();
      for(; cur!=end ; ++cur) {
	const StateDoublePair &sdp=*cur;
	STATE l=sdp.state;
	if(l==0) continue;
	double transP=sdp.logP;
	double sum=0.0;
	for(int i=0 ; i<L ; ++i)
	  sum+=exp(F(k,i)+transP+
		   hmm.getEmissionProb(l,S[i])+B(l,i+1)-seqP);
	A[k][l]+=seqWeight*sum;
      }
    }
    for(int k=1 ; k<Nq ; ++k)
      A[k][0]+=seqWeight*exp(F(k,L)+hmm.getLogTransProb(k,0)-seqP);
    A[0][0]=0;

    // Update discrete emission counts
    for(STATE q=1 ; q<Nq ; ++q) {
      for(int i=0 ; i<numDiscrete ; ++i) {
	const NmerChain &chain=bw.nmerChains[i];
	BOOM::Vector<double> V;
	const int firstPos=0;//### hmm.getOrder();
	for(int pos=firstPos ; pos<L ; ++pos)
	  V.push_back(F(q,pos+1)+B(q,pos+1));// ### added the +1's
	double den=sumLogProbs(V);
	for(int pos=firstPos ; pos<L ; ++pos) { 
	  const Emission &Si=S[pos];
	  /*const*/ double rho=F(q,pos+1)+B(q,pos+1)-den;// ### added the +1's
	  rho+=log(seqWeight);
	  NmerSymbol nmer=S[pos].getDiscrete(i);
	  Array1D<double> &nc=nmerCounts[q][i];
	  int ncSize=nc.size();
	  while(1) {
	    if(int(nmer)<ncSize) {
	      double &x=nc[nmer];
	      x=sumLogProbs(x,rho);
	    }
	    if(nmer==0) break;
	    nmer=chain.getSuffix(nmer);
	  }
	}
      }	
    }
  }//foreach sequence k
}



/****************************************************************
                          BaumThread2 methods
 ****************************************************************/
BaumThread2::BaumThread2(const BaumWelchMT &parent,
			 const BOOM::Vector<int> &seqs,
			 const int m,const int D,
			 const Array1D<float> &seqWeights)
  : bw(parent), seqs(seqs), sgn(D,D,m), m(m), D(D), 
    trainingSet(parent.trainingSet), hmm(parent.hmm),
    hmmGraph(parent.hmmGraph), edges(parent.edges), mu(parent.mu),
    seqWeights(seqWeights)
{
  sgn.setAllTo(0.0);
}



void BaumThread2::f()
{
  // Update covariance matrix
  const int Nq=bw.Nq;
  BOOM::Vector<int>::const_iterator cur=seqs.begin(), end=seqs.end();
  for(; cur!=end ; ++cur) {
    const int k=*cur;
    const EmissionSequence &S=*trainingSet[k];
    const int L=S.length();
    ForwardAlgorithm F(hmmGraph,S);
    BackwardAlgorithm B(hmmGraph,S);
    for(int j=0 ; j<m ; ++j) {
      for(STATE q=1 ; q<Nq ; ++q) {
	for(int i=0 ; i<L ; ++i) {
	  const Emission &Si=S[i];
	  const double rho=bw.computeRho(q,i,j,F,B,Si);
	  GSL::Vector SminusMu(D);
	  const GSL::Vector &mu_j=mu[j];
	  Si.getContinuous().sub(mu_j,SminusMu);
	  EdgeVector::const_iterator cur=edges.begin(), end=edges.end();
	  for(; cur!=end ; ++cur) {
	    const Edge &edge=*cur;
	    VertexId d=edge.first, c=edge.second;
	    const double Mdc=SminusMu[d]*SminusMu[c];
	    sgn[d][c][j]+=rho*Mdc;
	  }
	}
      }
    }//foreach component j...
  }//foreach sequence k
}



