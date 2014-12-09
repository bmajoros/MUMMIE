/*
 parse.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
*/
#include <math.h>
#include <iostream>
#include "BOOM/CommandLine.H"
#include "BOOM/Vector.H"
#include "BOOM/Sequence.H"
#include "BOOM/Array2D.H"
#include "BOOM/Regex.H"
#include "BOOM/File.H"
#include "BOOM/SummaryStats.H"
#include "BOOM/Constants.H"
#include "FastViterbi.H"
#include "EmissionLoader.H"
#include "IntervalDecoder.H"
using namespace std;


class Application
{
  Schema schema;
  int numContinuous;
  Regex deflineRegex, rangeRegex;
  String inHmmFile;
  bool wantPost, wantGff, justPosteriors;
  void run(const String &inSeqFile,HMM &,CommandLine &cmd);
  void emit(const String &substrate,int begin,int end,
	    double score,BOOM::Vector<double> *tracks,
	    EmissionSequence &seq,HMM &hmm,BOOM::Vector<int> &path);
  void addSeq(EmissionSequence &seq,GffFeature &f,HMM &hmm);
public:
  Application();
  int main(int argc,char *argv[]);
};


int main(int argc,char *argv[])
  {
    try
      {
	Application app;
	return app.main(argc,argv);
      }
    catch(const char *p)
      {
	cerr << p << endl;
      }
    catch(const string &msg)
      {
	cerr << msg.c_str() << endl;
      }
    catch(const exception &e)
      {
	cerr << "STL exception caught in main:\n" << e.what() << endl;
      }
    catch(...)
      {
	cerr << "Unknown exception caught in main" << endl;
      }
    return -1;
  }



Application::Application()
  : deflineRegex("\\s*>(\\S+)"), rangeRegex("(\\S+)-(\\S+)")
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"pg:PdI:c:");
  if(cmd.numArgs()!=2)
    throw string(
"parse [options] <*.hmm> <data.fastb>\n\
   where: -p = use posterior Viterbi\n\
          -g 2,3-7,8,13 = viterbi gff, foreground is defined by state list\n\
          -I 2,3-7,8,13 = score all contiguous runs of foreground states\n\
          -c 1e-10 = apply posterior cutoff when using -g or -I\n\
          -P = just emit posterior probabilities for all states\n\
          -d = second parameter is a directory rather than a file\n\
");
  inHmmFile=cmd.arg(0);
  String inSeqFile=cmd.arg(1);
  wantPost=cmd.option('p');
  wantGff=cmd.option('g');
  justPosteriors=cmd.option('P');
  if(wantGff && !wantPost) throw "-g requires -p";
  bool dir=cmd.option('d');

  // Load HMM and perform some initialization
  HMM hmm(inHmmFile);
  schema=hmm.getSchema();
  numContinuous=schema.getNumContinuous();

  // Load sequence(s)
  if(dir) {
    BOOM::Vector<String> files;
    File::getFileList(inSeqFile,files);
    BOOM::Vector<String>::iterator cur=files.begin(), end=files.end();
    for(; cur!=end ; ++cur) run(inSeqFile+'/'+*cur,hmm,cmd);
  }
  else run(inSeqFile,hmm,cmd);

  return 0;
}



void Application::run(const String &inSeqFile,HMM &hmm,CommandLine &cmd)
{
  Regex r("([^\/]+)\\\\.fastb$");
  if(!r.search(inSeqFile)) return;//throw inSeqFile;
  String substrate=r[1];
  EmissionLoader loader(hmm.getSchema());
  int numStates=hmm.countStates();
  EmissionSequence &seq=*loader.load(inSeqFile);
  seq.recode(hmm.getOrder()+1); // ###
  const int seqLen=seq.length();
  double cutoff=
    cmd.option('c') ? cmd.optParm('c').asFloat() : NEGATIVE_INFINITY;

  if(cmd.option('I')) {
    Set<int> &states=hmm.getForegroundStates();
    states.clear();
    BOOM::Vector<String> *fields=cmd.optParm('I').getFields(",");
    int n=fields->size();
    for(int i=0 ; i<n ; ++i) {
      String field=(*fields)[i];
      if(field.length()==0) continue;
      int begin, end;
      if(rangeRegex.match(field)) {
	begin=rangeRegex[1].asInt();
	end=rangeRegex[2].asInt();
      }
      else { begin=field.asInt(); end=begin; }
      for(int i=begin ; i<=end ; ++i) states.insert(i);
    }
    delete fields;
    IntervalDecoder D(hmm);
    BOOM::Vector<GffFeature> &features=*D.decode(seq,substrate,"MUMMIE",
						 "site");
    n=features.size();
    for(int i=0 ; i<n ; ++i) {
      GffFeature &f=features[i];
      if(f.getScore()>=cutoff) {
	addSeq(seq,f,hmm);
	/*
	EmissionSequence *sub=seq.getSubsequence(f.getBegin(),
						 f.getEnd()-f.getBegin());
	sub->unencode(hmm.getOrder());
	Sequence *subseq=sub->getDiscreteSeq(0);
	delete sub;
	String str;
	subseq->toString(hmm.getSchema().getAlphabet(0),0,subseq->getLength(),
			 str);
	delete subseq;
	f.getExtraFields().push_back(String("/seq=")+str);
	*/
	cout<<features[i];
      }
    }
    delete &features;
  }
  else if(justPosteriors) {
    HMMGraph hmmGraph(hmm);
    ForwardAlgorithm F(hmmGraph,seq);
    BackwardAlgorithm B(hmmGraph,seq);
    const double logP=B.getLogP();
    for(int pos=0 ; pos<seqLen ; ++pos)
      for(int q=1 ; q<numStates ; ++q)
	cout<<pos<<" "<<q<<" "<<exp(F(q,pos+1)+B(q,pos+1)-logP)<<endl;
  }
  else {
    FastViterbi viterbi(hmm,wantPost);
    Array1D<bool> foreground(numStates);
    if(wantGff) {
      foreground.setAllTo(false);
      BOOM::Vector<String> *fields=cmd.optParm('g').getFields(",");
      int n=fields->size();
      for(int i=0 ; i<n ; ++i) {
	String field=(*fields)[i];
	if(field.length()==0) continue;
	int begin, end;
	if(rangeRegex.match(field)) {
	  begin=rangeRegex[1].asInt();
	  end=rangeRegex[2].asInt();
	}
	else { begin=field.asInt(); end=begin; }
	for(int i=begin ; i<=end ; ++i) foreground[i]=true;
      }
      delete fields;
    }
    double pathScore;
    BOOM::Vector<double> stateScores;
    BOOM::Vector<int> &path=
      wantGff ? *viterbi.getPath(seq,pathScore,&stateScores)
      : *viterbi.getPath(seq,pathScore);
    cout<<"#pathScore="<<pathScore<<endl;
    int L=path.size();
    if(wantGff) {
      int begin;
      bool inSite=false;
      BOOM::Vector<double> scores;
      BOOM::Vector<double> tracks[numContinuous];
      for(int i=0 ; i<L ; ++i) {
	STATE q=path[i];
	double score=stateScores[i];
	if(foreground[q]) {
	  if(!inSite) { inSite=true; begin=i; }
	  scores.push_back(exp(score));
	  for(int j=0 ; j<numContinuous ; ++j)
	    tracks[j].push_back(seq[i][j]);
	}
	else {
	  if(inSite) { 
	    inSite=false;
	    SummaryStats stats(scores);
	    double score=stats.getMean();
	    scores.clear();
	    int end=i;
	    if(begin+1<end && score>=cutoff)
	      emit(substrate,begin,end,score,tracks,seq,hmm,path);
	    //cout<<substrate<<"\tbinding\tsite\t"<<begin+1<<"\t"<<end<<"\t"
	    //  <<score<<"\t+\t.\n";
	  }
	}
      }
      if(inSite) { 
	inSite=false;
	SummaryStats stats(scores);
	double score=stats.getMean();
	scores.clear();
	int end=L;
	if(begin+1<end && score>=cutoff)
	  emit(substrate,begin,end,score,tracks,seq,hmm,path);
	//cout<<substrate<<"\tbinding\tsite\t"<<begin+1<<"\t"<<end<<"\t"
	//    <<score<<"\t+\t.\n";
      }
    }
    else for(int i=0 ; i<L ; ++i) cout<<path[i]<<endl;
    delete &path;
  }
  delete &seq;
}



void Application::emit(const String &substrate,int begin,int end,
		       double score,BOOM::Vector<double> *tracks,
		       EmissionSequence &seq,HMM &hmm,
		       BOOM::Vector<int> &path)
{
  cout<<substrate<<"\tbinding\tsite\t"<<begin+1<<"\t"<<end<<"\t"
      <<score<<"\t+\t.\t";
  for(int i=0 ; i<numContinuous ; ++i) {
    SummaryStats stats(tracks[i]);
    double mean=stats.getMean();
    tracks[i].clear();
    cout<<schema.getContinuousName(i)<<"="<<mean<<";";
  }
  EmissionSequence *sub=seq.getSubsequence(begin,end-begin);
  sub->unencode(hmm.getOrder()+1);
  Sequence *subseq=sub->getDiscreteSeq(0);
  delete sub;
  String str;
  subseq->toString(hmm.getSchema().getAlphabet(0),0,subseq->getLength(),
		   str);
  delete subseq;
  cout<<"seq="<<str<<";path=";
  for(int pos=begin ; pos<end ; ++pos) {
    cout<<path[pos];
    if(pos+1<end) cout<<",";
  }
  cout<<";"<<endl;
}


void Application::addSeq(EmissionSequence &seq,GffFeature &f,HMM &hmm)
{
  EmissionSequence *sub=seq.getSubsequence(f.getBegin(),f.getLength());
  sub->unencode(hmm.getOrder()+1);
  Sequence *subseq=sub->getDiscreteSeq(0);
  delete sub;
  String str;
  subseq->toString(hmm.getSchema().getAlphabet(0),0,subseq->getLength(),
		   str);
  delete subseq;
  f.getExtraFields().push_back(String("seq=")+str);
}




