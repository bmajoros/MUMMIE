/****************************************************************
 SequenceSet.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "SequenceSet.H"
#include "BOOM/File.H"
#include "BOOM/SummaryStats.H"
#include "BOOM/GSL/GaussianDistribution.H"
#include "BOOM/DnaAlphabet.H"
#include "EmissionLoader.H"
using namespace std;
using namespace BOOM;


SequenceSet::SequenceSet()
{
  // ctor
}



int SequenceSet::getDim() const
{
  return (*seqs[0])[0].getContinuousDim();
}



void SequenceSet::normalize()
{
  const int dim=getDim(), N=size();
  transforms.resize(dim);
  for(int i=0 ; i<dim ; ++i) {
    BOOM::Vector<double> x;
    for(int j=0 ; j<N ; ++j) {
      EmissionSequence &S=*seqs[j];
      const int L=S.length();
      for(int k=0 ; k<L ; ++k) {
	x.push_back(S[k][i]);
      }
    }
    SummaryStats stats(x);
    double mean=stats.getMean(), var=stats.getVar(), m=stats.getMax();
    //const double slope=1/var;
    const double slope=100/m;
    const double intercept=-mean;
    transforms[i]=LinearFunc(slope,intercept,0);
  }	
  normalize(transforms);
}



void SequenceSet::normalize(const BOOM::Vector<LinearFunc> &t)
{
  const int dim=getDim(), N=size();
  for(int i=0 ; i<dim ; ++i) {
    const LinearFunc &f=transforms[i];
    for(int j=0 ; j<N ; ++j) {
      EmissionSequence &S=*seqs[j];
      const int L=S.length();
      for(int k=0 ; k<L ; ++k) {
	double &x=S[k].getContinuous()[i];
	x=f(x);
      }
    }
  }	
  if(&transforms!=&t) transforms=t;
}



void SequenceSet::saveTransforms(ofstream &os) const
{
  int n=transforms.size();
  os<<n<<endl;
  for(int i=0 ; i<n ; ++i) os<<transforms[i]<<endl;
}



void SequenceSet::loadTransforms(ifstream &is)
{
  int n;
  is>>n;
  transforms.resize(n);
  for(int i=0 ; i<n ; ++i) is>>transforms[i];
}



int SequenceSet::size() const
{
  return seqs.size();
}



void SequenceSet::resize(int s)
{
  seqs.resize(s);
}



EmissionSequence *SequenceSet::operator[](int i)
{
  return seqs[i];
}



const EmissionSequence *SequenceSet::operator[](int i) const
{
  return seqs[i];
}



void SequenceSet::load(const String &dir,Schema &schema,int maxFiles)
{
  BOOM::Vector<String> files, tmp;
  File::getFileList(dir,files);
  int numFiles=files.size();
  if(maxFiles>=0 && numFiles>maxFiles) numFiles=maxFiles;
  for(int i=0 ; i<numFiles ; ++i) {
    String filename=files[i];
    int L=filename.length();
    if(L<7) continue;
    String extension=filename.substr(L-5);
    extension.toupper();
    if(extension=="FASTB") tmp.push_back(filename);
  }
  files=tmp;
  numFiles=files.size();
  if(numFiles==0) throw "No files with extension \".fastb\" found";
  EmissionLoader loader(schema);
  for(int i=0 ; i<numFiles ; ++i) {
    const String fullpath=dir+"/"+files[i];
    EmissionSequence *seq=loader.load(fullpath);
    if(seq->length()==0) 
      throw String("ERROR: empty sequence in training set: ")+fullpath;
    seq->setFilename(fullpath);
    seqs.push_back(seq);
  }
}



void SequenceSet::save()
{
  const int n=seqs.size();
  for(int i=0 ; i<n ; ++i) {
    EmissionSequence *seq=seqs[i];
    seq->save(seq->getFilename());
  }
}



void SequenceSet::recode(int N)
{
  int numSeqs=seqs.size();
  for(int i=0 ; i<numSeqs ; ++i) seqs[i]->recode(N);
}



void SequenceSet::recode(int trackNum,int N)
{
  int n=seqs.size();
  for(int i=0 ; i<n ; ++i) seqs[i]->recode(trackNum,N);
}



void SequenceSet::addNoise(double mean,double var)
{
  GSL::GaussianDistribution gauss(mean,var);
  int N=seqs.size();
  for(int i=0 ; i<N ; ++i) {
    EmissionSequence &S=*seqs[i];
    const int L=S.length();
    for(int pos=0 ; pos<L ; ++pos) {
      Emission &e=S[pos];
      int D=e.getContinuousDim();
      GSL::Vector &v=e.getContinuous();
      for(int j=0 ; j<D ; ++j)
	v[j]+=gauss.random();
    }
  }
  
}


