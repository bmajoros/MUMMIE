/****************************************************************
 Kmeans.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "Kmeans.H"
#include "BOOM/Random.H"
#include "BOOM/Constants.H"
using namespace std;
using namespace BOOM;

Kmeans::Kmeans(int K,KmeansComparator &cmp)
  : K(K), cmp(cmp), centroids(K)
{
}



double Kmeans::cluster(const BOOM::Vector<GSL::Vector> &data,
		       int maxIter,double minChange,int numRestarts)
{
  ++numRestarts;
  double bestRSS;
  Array1D<GSL::Vector> bestCentroids;
  IntArray1D bestMembership;
  for(int i=0 ; i<numRestarts ; ++i) {
    double RSS=cluster(data,maxIter,minChange);
    if(i==0 || RSS<bestRSS) {
      bestRSS=RSS;
      bestCentroids=centroids;
      bestMembership=membership;
    }
  }
  centroids=bestCentroids;
  membership=bestMembership;
  iterationReport(data);
  return bestRSS;
}



double Kmeans::cluster_swap(const BOOM::Vector<GSL::Vector> &data,
			    int maxIter,double minChange,int numSwaps)
{
  const int N=data.size();
  membership.resize(N);
  const int dim=data[0].getDim();
  allocateCentroids(dim);
  cout<<"initializing clusters..."<<endl;
  randomAssignment(data);
  ++numSwaps;
  double bestRSS;
  Array1D<GSL::Vector> bestCentroids;
  IntArray1D bestMembership;
  for(int i=0 ; i<numSwaps ; ++i) {
    double RSS=cluster_core(data,maxIter,minChange);
    if(i==0 || RSS<bestRSS) {
      bestRSS=RSS;
      bestCentroids=centroids;
      bestMembership=membership;
    }
    perturb();
  }
  centroids=bestCentroids;
  membership=bestMembership;
  iterationReport(data);
  return bestRSS;
}



void Kmeans::perturb()
{
  const int N=membership.size();
  const int victim=RandomNumber(N);
  membership[victim]=(membership[victim]+RandomNumber(K-1)+1)%K;
}



double Kmeans::cluster(const BOOM::Vector<GSL::Vector> &data,
		       int maxIter,double minChange)
{
  const int N=data.size();
  membership.resize(N);
  const int dim=data[0].getDim();
  allocateCentroids(dim);
  cout<<"initializing clusters..."<<endl;
  randomAssignment(data);
  cluster_core(data,maxIter,minChange);
}



double Kmeans::cluster_core(const BOOM::Vector<GSL::Vector> &data,
			    int maxIter,double minChange)
{
  bool changes=true;
  int iter=1;
  double prevRSS=-1;
  while(changes) {
    cout<<"ITERATION #"<<iter;
    Array1D<GSL::Vector> debug=centroids;
    updateCentroids(data);
    const double RSS=computeRSS(data);
    cout<<"\tRSS="<<RSS<<endl;
    iterationReport(data);
    changes=reassign(data);
    ++iter;
    if(prevRSS>0) {
      const double deltaRSS=fabs(RSS-prevRSS);
      if(deltaRSS<=minChange) {	
	prevRSS=RSS; 
	break; 
      }
    }
    prevRSS=RSS;
    if(maxIter>0 && iter>=maxIter) {
      break;
    }
  }
  return prevRSS;
}



void Kmeans::iterationReport(const BOOM::Vector<GSL::Vector> &data)
{
  IntArray1D counts(K);
  counts.setAllTo(0);
  const int N=membership.size();
  for(int i=0 ; i<N ; ++i) ++counts[membership[i]];
  cout<<"cluster sizes: ";
  for(int i=0 ; i<K ; ++i) cout<<counts[i]<<" ";
  cout<<endl;
  return;

  cout.precision(3);
  cout<<"centroids:\n";
  const int dim=centroids[0].getDim();
  for(int i=0 ; i<K ; ++i) {
    cout<<"\t";
    printVector(centroids[i]);
    cout<<endl;
    IntVector members;
    getClusterMembers(i,members);
    for(int i=0 ; i<members.size() ; ++i) {
      cout<<"\t";
      printVector(data[members[i]]);
      cout<<endl;
    }
    cout<<endl;
  }
}



void Kmeans::printVector(const GSL::Vector &v)
{
  const int dim=v.getDim();
  for(int j=0 ; j<dim ; ++j) {
    cout<<v[j];
    if(j<dim-1) cout<<"\t";
  }
}



void Kmeans::getClusterMembers(int clusterIndex,IntVector &into)
{
  into.clear();
  const int N=membership.size();
  for(int i=0 ; i<N ; ++i)
    if(membership[i]==clusterIndex)
      into.push_back(i);
}



void Kmeans::allocateCentroids(int dim)
{
  for(int i=0 ; i<K ; ++i) centroids[i].resize(dim);
}



bool Kmeans::reassign(const BOOM::Vector<GSL::Vector> &data)
{
  bool changes=false;
  const int N=data.size();
  for(int i=0 ; i<N ; ++i) {
    const GSL::Vector &record=data[i];
    int bestIndex=membership[i];
    double bestDistance=cmp(record,centroids[bestIndex]);
    for(int j=0 ; j<K ; ++j) {
      double distance=cmp(record,centroids[j]);
      if(isFinite(distance) && 
	 (distance<bestDistance ||
	  distance==bestDistance && j<bestIndex)) {
	bestIndex=j;
	bestDistance=distance;
      }
    }
    if(bestIndex!=membership[i]) {
      membership[i]=bestIndex;
      changes=true;
    }
  }
  return changes;
}



void Kmeans::updateCentroids(const BOOM::Vector<GSL::Vector> &data)
{
  IntArray1D counts(K); counts.setAllTo(0);
  for(int i=0 ; i<K ; ++i) centroids[i].setAllTo(0.0);
  const int N=data.size();
  for(int i=0 ; i<N ; ++i) {
    const int memb=membership[i];
    centroids[memb].accumulate(data[i]);
    ++counts[memb];
  }
  const int dim=data[i].getDim();
  for(int i=0 ; i<K ; ++i) {
    GSL::Vector &centroid=centroids[i];
    const int size=counts[i];
    for(int j=0 ; j<dim ; ++j) centroid[j]/=size;
  }
}



void Kmeans::randomAssignment(const BOOM::Vector<GSL::Vector> &data)
{
  const int N=data.size();
  for(int i=0 ; i<N ; ++i) membership[i]=RandomNumber(K);
}



const Array1D<GSL::Vector> &Kmeans::getCentroids() const
{
  return centroids;
}



const IntArray1D &Kmeans::getMembership() const
{
  return membership;
}



double Kmeans::computeRSS(const PointSet &data)
{
  double SS=0;
  const int N=data.size();
  for(int i=0 ; i<N ; ++i) {
    const GSL::Vector x=data[i];
    const GSL::Vector centroid=centroids[membership[i]];
    const double dist=cmp(x,centroid);
    if(!isFinite(dist)) {cout<<"INFINITE DISTANCE: "<<x<<" VS "<<centroid<<endl;INTERNAL_ERROR;}
    SS+=dist*dist;
  }
  const double RSS=sqrt(SS);
  return RSS;
}



double KmeansEuclidean::operator()(const GSL::Vector &a,const GSL::Vector &b)
{
  const int D=a.getDim();
  double dist=0.0;
  for(int i=0 ; i<D ; ++i) {
    const double delta=a[i]-b[i];
    dist+=delta*delta;
  }
  return sqrt(dist);
}



double KmeansAnticorrelation::operator()(const GSL::Vector &a,
					 const GSL::Vector &b)
{
  throw "KmeansAnticorrelation::operator() not implemented";

  // 1/|r|
}



