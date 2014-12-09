/*
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 */
#include <string>
#include <iostream>
#include <fstream>
#include "BOOM/Random.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Vector.H"
#include "BOOM/SparseGraph.H"
#include "BOOM/Constants.H"
#include "BOOM/SummaryStats.H"
#include "BOOM/Regex.H"
#include "BaumWelch.H"
#include "BaumWelchMT.H"
#include "EmissionLoader.H"
#include "SequenceSet.H"
#include "TieProfile.H"
using namespace std;

const int WINDOW_SIZE=100;

class Application
{
  Alphabet alphabet;
  int sigID, openID;
  SequenceSet trainingSet;
  bool findMotif(int peakPos,int &from,int &to,
		 const BOOM::Vector<String> &motifs,int margin,
		 EmissionSequence &S);
  bool match(const String &pattern,EmissionSequence &S,int begin,int len);
public:
  int go(int argc,char *argv[]);
};


int main(int argc,char *argv[])
{
  try
    {
      Application app;
      return app.go(argc,argv);
    }
  catch(const char *msg)
    {
      cerr << msg << endl;
      return -1;
    }
  catch(const String &s)
    {
      cerr << s.c_str() << endl;
      return -1;
    }
  catch(...)
    {
      cerr << "unknown exception caught in main()" << endl;
      return -1;
    }
  return 0;
}



int Application::go(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=5) 
    throw "\n\
hack-known-sites [options] <training-dir> <schema.txt> <signal-name> <margin-around-peak> <motifs>\n\
";
  String dir=cmd.arg(0);
  String schemaFile=cmd.arg(1);
  String signalName=cmd.arg(2);
  int margin=cmd.arg(3).asInt();
  String motifs=cmd.arg(4);

  // Misc. initialization
  cout.precision(10);
  cerr.precision(10);
  BOOM::catchFloatOverflow();
  BOOM::Vector<String> &motifV=*motifs.getFields(",");
  Regex r("([^\/]+)\\.fastb");

  // Load schema
  Schema schema(schemaFile);
  sigID=schema.lookupContinuousID(signalName);
  openID=schema.lookupContinuousID("unpaired");
  alphabet=schema.getAlphabet(0);

  // Load the training set
  cerr<<"loading training data"<<endl;
  trainingSet.load(dir,schema);

  // Process sequences
  cerr<<"processing sequences"<<endl;
  BOOM::Vector<float> all, cluster;
  BOOM::Array1D< BOOM::Vector<double> > offsets(2*WINDOW_SIZE+1);
  BOOM::Vector<double> motifClusterSizes, otherClusterSizes;
  const int N=trainingSet.size();
  int nextSeqID=1;
  for(int i=0 ; i<N ; ++i) {
    EmissionSequence &S=*trainingSet[i];
    if(!r.search(S.getFilename())) throw S.getFilename();
    String substrate=r[1];
    int begin;
    double prevX=0;
    double maxX=-1;
    int maxPos;
    const int L=S.length();
    for(int pos=0 ; pos<L ; ++pos) {
      const Emission &e=S[pos];
      const double x=e[sigID];
      all.push_back(e[openID]);
      if(prevX==0 && x>0) { begin=pos; maxX=x; maxPos=pos; }
      else if(prevX>0 && x>0 && x>maxX) { maxX=x; maxPos=pos; }
      else if(prevX>0 && x==0) {
	int from, to;
	if(1) {
	  if(findMotif(maxPos,from,to,motifV,margin,S)) {
	  /*
	  cout<<substrate<<"\tpositive\tgroup\t"<<begin<<"\t"<<pos<<"\t1\t+\t.\n";
	  cout<<substrate<<"\tmotif\tsite\t"<<from<<"\t"<<to<<"\t1\t+\t.\n";
	  */
	  //cout<<maxX<<endl;
	  /*
	  EmissionSequence *subseq=S.getSubsequence(begin,pos-begin);
	  Sequence *seq=subseq->getDiscreteSeq(0);
	  delete subseq;
	  cout<<">"<<nextSeqID++<<endl;
	  seq->printOn(cout,schema.getAlphabet(0));
	  cout<<endl;
	  delete seq;
	  */

	  for(int delta=-WINDOW_SIZE ; delta<=WINDOW_SIZE ; ++delta) {
	    //int x=to+delta;
	    int x=maxPos+delta;
	    if(x<0 || x>=L) continue;
	    double y=S[x][openID];
	    offsets[delta+WINDOW_SIZE].push_back(y);
	  }

	  //motifClusterSizes.push_back(pos-begin);
	  }
	}
	else {//for(int i=begin ; i<pos ; ++i) {
	  /*
	  cout<<substrate<<"\tnegative\tgroup\t"<<begin+1<<"\t"<<pos<<"\t1\t+\t.\n";
	  cout<<substrate<<"\tdecoy\tsite\t"<<maxPos-margin<<"\t"<<maxPos+margin<<"\t1\t+\t.\n";
	  */
	  //cout<<maxX<<endl;
	  /*
	  //EmissionSequence *subseq=S.getSubsequence(begin,pos-begin);
	  if(maxPos>margin && maxPos+margin<L) {
	    EmissionSequence *subseq=S.getSubsequence(maxPos-margin,2*margin);
	    Sequence *seq=subseq->getDiscreteSeq(0);
	    delete subseq;
	    cout<<">"<<nextSeqID++<<endl;
	    seq->printOn(cout,schema.getAlphabet(0));
	    cout<<endl;
	    delete seq;
	  }
	  */
	  for(int delta=-WINDOW_SIZE ; delta<=WINDOW_SIZE ; ++delta) {
	    int x=maxPos+delta;
	    if(x<0 || x>=L) continue;
	    double y=S[x][openID];
	    offsets[delta+WINDOW_SIZE].push_back(y);
	  }

	  //cluster.push_back(S[i][openID]);
	  int start=maxPos-margin, stop=maxPos+margin+1;
	  if(start<0) start=0;
	  if(stop>L) stop=L;
	  else for(int i=start ; i<stop ; ++i) 
	    cluster.push_back(S[i][openID]);
	  otherClusterSizes.push_back(pos-begin);
	}
      }
      prevX=x;
    }
  }
  /*
  SummaryStats ss(all), ssC(cluster);
  SummaryStats sizes1(motifClusterSizes), sizes2(otherClusterSizes);
  cerr<<"overall stats: "<<ss.getMean()<<" +/- "<<ss.getStdDev()<<endl;
  cerr<<"cluster stats: "<<ssC.getMean()<<" +/- "<<ssC.getStdDev()<<" N="<<ssC.getN()<<endl;
  cerr<<"motif cluster sizes: "<<sizes1.getMean()<<" +/- "<<sizes1.getStdDev()<<" N="<<sizes1.getN()<<endl;
  cerr<<"non-motif cluster sizes: "<<sizes2.getMean()<<" +/- "<<sizes2.getStdDev()<<" N="<<sizes2.getN()<<endl;
  */
  for(int delta=-WINDOW_SIZE ; delta<=WINDOW_SIZE ; ++delta) {
    int x=delta+WINDOW_SIZE;
    SummaryStats ss(offsets[x]);
    cout<<delta<<"\t"<<ss.getMean()<<endl;
  }

  return 0;
}


bool Application::findMotif(int peakPos,int &from,int &to,
			    const BOOM::Vector<String> &motifs,int margin,
			    EmissionSequence &S) 
{
  int L=S.length();
  const int numMotifs=motifs.size();
  const int motifLen=motifs[0].length();
  int begin=peakPos-margin, end=peakPos+margin-motifLen+1;
  if(begin<0) begin=0;
  if(end>=L) end=L-1;
  for(int i=begin ; i<end ; ++i) {
    for(int j=0 ; j<numMotifs ; ++j)
      if(match(motifs[j],S,i,motifLen)) {
	from=i;
	to=i+motifLen;
	return true;
      }
  }
  return false;
}



bool Application::match(const String &pattern,EmissionSequence &S,
			int begin,int len)
{
  int end=begin+len;
  for(int i=begin ; i<end ; ++i)
    if(S[i].getDiscrete(0)!=alphabet.lookup(pattern[i-begin])) return false;
  return true;
}



