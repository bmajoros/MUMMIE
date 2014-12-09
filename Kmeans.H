/****************************************************************
 Kmeans.H
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_Kmeans_H
#define INCL_Kmeans_H
#include <iostream>
#include "BOOM/GSL/Vector.H"
#include "BOOM/Vector.H"
#include "BOOM/Array1D.H"
using namespace std;
using namespace BOOM;

struct KmeansComparator {
  virtual double operator()(const GSL::Vector &,const GSL::Vector &)=0;
};

struct KmeansEuclidean : public KmeansComparator {
  virtual double operator()(const GSL::Vector &,const GSL::Vector &);
};

struct KmeansAnticorrelation : public KmeansComparator {
  virtual double operator()(const GSL::Vector &,const GSL::Vector &);
};

class Kmeans {
public:
  typedef const BOOM::Vector<GSL::Vector> PointSet;
  Kmeans(int K,KmeansComparator &);
  double cluster(const PointSet &data,int maxIter,double minChange,
		 int numRestarts);
  double cluster_swap(const BOOM::Vector<GSL::Vector> &data,
		      int maxIter,double minChange,int numSwaps);
  const Array1D<GSL::Vector> &getCentroids() const;
  const IntArray1D &getMembership() const;
protected:
  const int K;
  KmeansComparator &cmp;
  Array1D<GSL::Vector> centroids;
  IntArray1D membership;
  void allocateCentroids(int dim);
  void randomAssignment(const PointSet &);
  void updateCentroids(const PointSet &);
  bool reassign(const PointSet &);
  void iterationReport(const PointSet &);
  void getClusterMembers(int clusterIndex,IntVector &);
  void printVector(const GSL::Vector &);
  double computeRSS(const PointSet &); // RSS measures clustering quality
  double cluster(const BOOM::Vector<GSL::Vector> &data,int maxIter,
		 double minChange);
  double cluster_core(const BOOM::Vector<GSL::Vector> &data,
		      int maxIter,double minChange);
  void perturb();
};

#endif

