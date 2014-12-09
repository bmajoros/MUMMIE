/**************************************************************
 classify.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/

#include <string>
#include <iostream>
#include <fstream>
#include "BOOM/FastaReader.H"
#include "BOOM/Sequence.H"
#include "BOOM/DnaAlphabet.H"
#include "BOOM/CommandLine.H"
#include "HMM.H"
using namespace std;

class Application
{
  DnaAlphabet alphabet;
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



void usage()
{
    throw "\n\
classify [options] <0.hmm> <1.hmm> <*.fasta>\n  \
    where [options] are:\n\
        -t <threshold> : use threshold for log-llh-ratio (default=0)\n\n";
}


int Application::go(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"t:");
  if(cmd.numArgs()!=3) usage();
  String hmmFile0=cmd.arg(0);
  String hmmFile1=cmd.arg(1);
  String fastaFile=cmd.arg(2);
  double threshold=0;
  if(cmd.option('t')) threshold=cmd.optParm('t').asDouble();
  
  // Load HMMs
  HMM hmm0(hmmFile0,alphabet), hmm1(hmmFile1,alphabet);

  // Process the FASTA file
  FastaReader reader(fastaFile);
  String defline, sequence;

  while(reader.nextSequence(defline,sequence))
    {
      Sequence seq(sequence,alphabet);
      double score0=hmm0.getLogP(seq);
      double score1=hmm1.getLogP(seq);
      double ratio=score1-score0;
      if(ratio>threshold) cout<<1<<endl;
      else cout<<0<<endl;
    }

  /*
  while(reader.nextSequence(defline,sequence))
    {
      Sequence seq(sequence,alphabet);
      Sequence s0, s1, s2;
      s0=seq;
      int L=seq.getLength();
      seq.getSubsequence(1,L-1,s1);
      seq.getSubsequence(2,L-2,s2);
      
      double s00=hmm0.getLogP(s0);
      double s01=hmm0.getLogP(s1);
      double s02=hmm0.getLogP(s2);
      double score0=s00; if(s01>s00) score0=s01; if(s02>score0) score0=s02;
      
      double s10=hmm1.getLogP(s0);
      double s11=hmm1.getLogP(s1);
      double s12=hmm1.getLogP(s2);
      double score1=s10; if(s11>s10) score1=s11; if(s12>score1) score1=s12;
      
      double ratio=score1-score0;
      if(ratio>threshold) cout<<1<<endl;
      else cout<<0<<endl;
    }
  */
  
  return 0;
}


