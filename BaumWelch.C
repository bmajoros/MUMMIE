


OBSOLETE.




/**************************************************************
BaumWelch.C : the Baum-Welch Expectation Maximization algorithm
bmajoros@duke.edu

This is OPEN SOURCE SOFTWARE governed by the ARTISTIC LICENSE.
***************************************************************/
#include <fstream>
#include "BaumWelch.H"
#include "BOOM/SummaryStats.H"
#include "BOOM/Correlation.H"
#include "BOOM/GSL/GaussianDistribution.H"
#include "BOOM/Array1D.H"
#include "BOOM/Array2D.H"
#include "BOOM/Array3D.H"
#include "BOOM/Constants.H"
#include "BOOM/SumLogProbs.H"


/****************************************************************
                        BaumWelchMT methods
 ****************************************************************/

BaumWelch::BaumWelch(HMM &hmm,long maxIterations,
		     const BOOM::Vector<EmissionSequence*> &trainingSet,
		     SparseGraph *dependencyGraph,int maxSampleSize)
  : hmm(hmm), numHmmStates(hmm.countStates()),
    numTrain(trainingSet.size()), likelihoods(trainingSet.size()),
    trainingSet(trainingSet),  maxIterations(maxIterations),
    hmmGraph(hmm), schema(hmm.getSchema()), maxSampleSize(maxSampleSize),
    dependencyGraph(dependencyGraph)
{
  if(dependencyGraph) dependencyGraph->getEdges(edges);
  else throw "BaumWelch requires dependency graph";
  mainAlgorithm();
}



void BaumWelch::getGlobalStats()
{
  const int D=schema.getNumContinuous();
  globalMeans.resize(D);
  globalCov.resize(D,D);
  const int numTrain=trainingSet.size();

  // Compute means and variances
  for(int d=0 ; d<D ; ++d) {
    DoubleVector V;
    for(int i=0 ; i<numTrain ; ++i) {
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
    globalCov(d,d)=stats.getVar();
  }
  
  // Compute covariances
  if(dependencyGraph) {sparseStats();return;}
  for(int d1=0 ; d1<D ; ++d1) {
    double sd1=sqrt(globalCov(d1,d1));
    for(int d2=d1+1 ; d2<D ; ++d2) {
      //if(d2==d1) continue;
      double sd2=sqrt(globalCov(d2,d2));
      DoubleVector V1, V2;
      for(int i=0 ; i<numTrain ; ++i) {
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
      globalCov(d1,d2)=globalCov(d2,d1)=cor.getR()*sd1*sd2;
    }
  }
}



void BaumWelch::sparseStats()
{
  // Compute covariances
  const int numTrain=trainingSet.size();
  BOOM::Vector< pair<VertexId,VertexId> >::iterator cur=edges.begin(), 
    end=edges.end();
  for(; cur!=end ; ++cur) {
    pair<VertexId,VertexId> &p=*cur;
    VertexId d1=p.first, d2=p.second;
    if(d2==d1) continue;
    double sd1=sqrt(globalCov(d1,d1));
    double sd2=sqrt(globalCov(d2,d2));
    DoubleVector V1, V2;
    for(int i=0 ; i<numTrain ; ++i) {
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
    globalCov(d2,d1)=globalCov(d1,d2)=cor.getR()*sd1*sd2;
  }  
}


double BaumWelch::computeRho(STATE q,int i,int j,ForwardAlgorithm &F,
			     BackwardAlgorithm &B,const Emission &e)
{
  ++i;
  const int m=hmm.numMixtureComponents();
  GaussianMixture &mix=hmm.getEmissionDistr(q);
  const GSL::Vector &Si=e.getContinuous();
  BOOM::Vector<double> V;
  for(STATE s=1 ; s<numHmmStates ; ++s) 
    V.push_back(safeAdd(F(s,i),B(s,i)));
  double den=sumLogProbs(V);
  V.clear();
  for(int k=0 ; k<m ; ++k) 
    V.push_back(safeAdd(mix.getLogCoef(k),mix.getDistr(k).logDensity(Si)));
  den+=sumLogProbs(V);
  double num=
    F(q,i)+B(q,i)+mix.getLogCoef(j)+mix.getDistr(j).logDensity(Si);
  double rho=exp(num-den);
  if(!isFinite(rho)) {
    cout<<"m="<<m<<" Si="<<Si<<" i="<<i<<" j="<<j<<" q="<<q<<" den="<<den<<" num="<<num<<" F="<<F(q,i+1)<<" B="<<B(q,i+1)<<" lambda="<<mix.getCoef(j)<<" density="<<mix.getDistr(j).logDensity(Si)<<endl;
    for(STATE s=0 ; s<numHmmStates ; ++s) 
      cout<<"F("<<s<<",i+1)="<<F(s,i+1)<<" B("<<s<<",i+1)="<<B(s,i+1)<<endl;
    for(int k=0 ; k<m ; ++k) 
      cout<<"logCoef="<<log(mix.getCoef(k))<<" logDensity="
	  <<mix.getDistr(k).logDensity(Si)<<endl;
    INTERNAL_ERROR;
  }
  //cout<<"rho="<<rho<<" num="<<num<<" den="<<den<<endl;
  return rho;
}



void BaumWelch::mainAlgorithm()
{
  // Initialize some constants
  const int Nq=numHmmStates;
  const int m=hmm.numMixtureComponents();
  const int D=schema.getNumContinuous();
  const int K=trainingSet.size();

  // Observe global stats
  cout<<"initializing global stats"<<endl;
  getGlobalStats();
  Array1D<GaussianDistribution> gauss(D);
  for(int d=0 ; d<D ; ++d) 
    gauss[d]=GaussianDistribution(globalMeans[d],globalCov(d,d));

  // Initialize the means and covariance matrices
  Array1D<GSL::Vector> mu(m);
  for(int j=0 ; j<m ; ++j) {
    GSL::Vector &mu_j=mu[j];
    mu_j.resize(D);
    for(int d=0 ; d<D ; ++d) mu_j[d]=gauss[d].random();
  }
  for(int j=0 ; j<m ; ++j) {
    MultiGauss gauss(mu[j],globalCov);
    for(STATE q=1 ; q<Nq ; ++q) {
      GaussianMixture &mix=hmm.getEmissionDistr(q);
      mix.setDistr(j,gauss);
    }
  }

  // Initialize lambda
  for(STATE q=1 ; q<Nq ; ++q) {
    GaussianMixture &mix=hmm.getEmissionDistr(q);
    Array1D<double> lambda(m);
    double sum=0;
    for(int j=0 ; j<m ; ++j) sum+=lambda[j]=Random0to1();
    for(int j=0 ; j<m ; ++j) mix.setCoef(j,lambda[j]/sum);
  }

  // Initialize transition counts
  A.resize(Nq,Nq);

  //====================================================== EM EM EM EM
  // Perform EM
  //ofstream osRho("rho.log");//###
  cout<<"performing EM"<<endl;
  DblArray1D msd(m);
  DblArray2D mun(D,m), ln(Nq,m), ld(Nq,m);
  DblArray3D sgn(D,D,m);
  double prevLL;
  for(int iter=0 ; iter<maxIterations ; ++iter) {
    double LL=0;
    cout<<"ITERATION #"<<iter<<"\t";cout.flush();
    //osRho<<"ITERATION #"<<iter<<endl;
    msd.setAllTo(0.0);mun.setAllTo(0.0);ln.setAllTo(0.0);
    ld.setAllTo(0.0);A.setAllTo(0);
    hmmGraph.updateProbs();
    for(int k=0 ; k<K ; ++k) {
      const EmissionSequence &S=*trainingSet[k];
      const int L=S.length();
      ForwardAlgorithm F(hmmGraph,S);
      BackwardAlgorithm B(hmmGraph,S);
      for(int j=0 ; j<m ; ++j) {
	for(STATE q=1 ; q<Nq ; ++q) {
	  for(int i=0 ; i<L ; ++i) {
	    const Emission &Si=S[i];
	    const double rho=computeRho(q,i,j,F,B,Si);
	    //osRho<<rho<<endl;
	    //if(!isFinite(rho) || rho==0) {cout<<"rho="<<rho<<" Si="<<Si.getContinuous()<<endl;INTERNAL_ERROR;}
	    ln[q][j]+=rho;
	    msd[j]+=rho;
	    for(int k=0 ; k<m ; ++k) ld[q][j]+=computeRho(q,i,k,F,B,Si);
	    for(int d=0 ; d<D ; ++d) mun[d][j]+=rho*Si[d];
	  }
	}
      }//foreach component j...

      // Update transition counts
      double seqP=F.getLogP();
      LL+=seqP;
      for(int k=0 ; k<Nq ; ++k) {
        for(int l=1 ; l<Nq ; ++l) {   
	  double sum=0.0;
	  for(int i=0 ; i<L ; ++i) {
	    sum+=exp(F(k,i)+log(hmm.getTransitionProb(k,l))+ 
		     hmm.getEmissionProb(l,S[i])+B(l,i+1)-seqP);
	    if(!isFinite(sum)) {cout<<"sum="<<sum<<endl;INTERNAL_ERROR;}
	  }
	  A[k][l]+=sum;
        }
      }
      for(int k=1 ; k<Nq ; ++k)
	A[k][0]+=exp(F(k,L)+log(hmm.getTransitionProb(k,0))-seqP);
      A[0][0]=0;
    }//foreach sequence k

    // Update transition probabilities
    for(int k=0 ; k<Nq ; ++k) {
      double sum=0;
      for(int l=0 ; l<Nq ; ++l) sum+=A[k][l];
      if(sum==0) for(int l=0 ; l<Nq ; ++l) 
	hmm.setTransitionProb(k,l,BOOM::NEGATIVE_INFINITY);
      else for(int l=0 ; l<Nq ; ++l) {
	double p=A[k][l]/sum;
	hmm.setTransitionProb(k,l,p);
      }
    }

    // Update mixture coefficients and means
    for(int j=0 ; j<m ; ++j) {
      for(STATE q=1 ; q<Nq ; ++q) {
	double lambda=ln[q][j]/ld[q][j];
	if(!isFinite(lambda)) lambda=0;
	hmm.getEmissionDistr(q).setCoef(j,lambda);
      }
      for(int d=0 ; d<D ; ++d) {
	double mu_j=mun[d][j]/msd[j];
	if(!isFinite(mu_j)) mu_j=0;
	mu[j][d]=mu_j;
	//if(!isFinite(mu[j][d])) {cout<<"mun="<<mun[d][j]<<" msd="<<msd[j]<<endl;;INTERNAL_ERROR;}
      }
    }

    // Update covariance matrix
    sgn.setAllTo(0.0);
    for(int k=0 ; k<K ; ++k) {
      const EmissionSequence &S=*trainingSet[k];
      const int L=S.length();
      ForwardAlgorithm F(hmmGraph,S);
      BackwardAlgorithm B(hmmGraph,S);
      for(int j=0 ; j<m ; ++j) {
	for(STATE q=1 ; q<Nq ; ++q) {
	  for(int i=0 ; i<L ; ++i) {
	    const Emission &Si=S[i];
	    const double rho=computeRho(q,i,j,F,B,Si);
	    if(!isFinite(rho)) {cout<<"rho="<<rho<<endl;INTERNAL_ERROR;}
	    GSL::Vector SminusMu(D);
	    const GSL::Vector &mu_j=mu[j];
	    Si.getContinuous().sub(mu_j,SminusMu);
	    EdgeVector::iterator cur=edges.begin(), end=edges.end();
	    for(; cur!=end ; ++cur) {
	      const Edge &edge=*cur;
	      VertexId d=edge.first, c=edge.second;
	      const double Mdc=SminusMu[d]*SminusMu[c];
	      //if(!isFinite(Mdc)) {cout<<"SminusMu="<<SminusMu<<" mu_j="<<mu_j<<endl;INTERNAL_ERROR;}
	      sgn[d][c][j]+=rho*Mdc;
	      if(!isFinite(sgn[d][c][j])) sgn[d][d][j]=0;
	      if(c!=d) sgn[c][d][j]=sgn[d][c][j];
	    }
	  }
	}
      }//foreach component j...
    }//foreach sequence k
    for(int j=0 ; j<m ; ++j) {
      GSL::Matrix Cj(D,D);
      Cj.setAllTo(0.0);
      EdgeVector::iterator cur=edges.begin(), end=edges.end();
      for(; cur!=end ; ++cur) {
	const Edge &edge=*cur;
	VertexId d=edge.first, c=edge.second;
	Cj(c,d)=Cj(d,c)=sgn[d][c][j]/msd[j];
      }
      if(validMean(mu[j]) && validCov(Cj)) {
	MultiGauss gauss(mu[j],Cj);
	//MultiGauss gauss(mu[j],globalCov);
	for(STATE q=1 ; q<Nq ; ++q)
	  hmm.getEmissionDistr(q).setDistr(j,gauss);
      }
      else {
	MultiGauss gauss;
	for(STATE q=1 ; q<Nq ; ++q) {
	  GaussianMixture &mix=hmm.getEmissionDistr(q);
	  mix.setDistr(j,gauss);
	  mix.setCoef(j,0);
	}
      }
    }

    cout<<"LL="<<LL;
    if(iter>0 && LL<prevLL) {
      cout<<"  <<<<< LL DECREASE";
      //osRho.close();
      //INTERNAL_ERROR;
    }
    cout<<endl;
    prevLL=LL;
  }
}



bool BaumWelch::validMean(const GSL::Vector &mu) 
{
  const int n=mu.getDim();
  for(int i=0 ; i<n ; ++i) if(!isFinite(mu[i])) return false;
  return true;
}



bool BaumWelch::validCov(const GSL::Matrix &M)
{
  const int n=M.getNumRows();
  for(int i=0 ; i<n ; ++i)
    for(int j=0 ; j<n ; ++j)
      if(!isFinite(M(i,j))) return false;
  return true;
}



