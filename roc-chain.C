/*
 roc-chain.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
*/
#include <math.h>
#include <iostream>
#include "BOOM/CommandLine.H"
#include "BOOM/Vector.H"#include "BOOM/Sequence.H"
#include "BOOM/Array2D.H"
#include "BOOM/Regex.H"
#include "BOOM/File.H"
#include "BOOM/SummaryStats.H"
#include "BOOM/Constants.H"
#include "BOOM/VectorSorter.H"
#include "BOOM/PureDnaAlphabet.H"
#include "EGGS/WMM.H"
#include "EGGS/GarbageCollector.H"
#include "FastViterbi.H"
#include "EmissionLoader.H"
#include "SequenceSet.H"
using namespace std;

BOOM::Alphabet alphabet=PureDnaAlphabet::global();

const int SIGNAL_SLICES=100; //50;
const int LLR_SLICES=100; //10;

struct Hit {
  float score;
  int begin, end;
  bool hit, miss, optimal;
  void reset() {hit=miss=optimal=false;}
  Hit(float s,int b,int e) : score(s), begin(b), end(e) {}
  bool overlap(const Hit &h) {return begin<h.end && end>h.begin;}
};

struct RocPoint {
  int TP, FP, TN, FN;
  float Sn, FPR, F, threshold;
  RocPoint() : TP(0), FP(0), TN(0), FN(0), Sn(NEGATIVE_INFINITY),
	       FPR(NEGATIVE_INFINITY), F(NEGATIVE_INFINITY) {}
  void init(int tp,int fp,int fn,int tn,float t=0) {
    TP=tp;    FP=fp;    FN=fn;    TN=tn;    threshold=t;
    Sn=TP+FN>0 ? TP/float(TP+FN) : 0.0;
    //Sp=TP+FP>0 ? TP/float(TP+FP) : 0.0;
    //FDR=FP+TP>0 ? FP/float(TP+FP) : 0.0;
    FPR=FP+TN>0 ? FP/float(FP+TN) : 0.0;
    F=2*TP/float(2*TP+FP+FN);
    //cout<<F<<" "<<Sn<<" "<<FPR<<" "<<threshold<<endl;
  }
  };

class RocComparator : public Comparator<RocPoint> {
public:
  virtual bool equal(RocPoint &a,RocPoint &b) {return a.FPR==b.FPR;}
  virtual bool greater(RocPoint &a,RocPoint &b) {return a.FPR>b.FPR;}
  virtual bool less(RocPoint &a,RocPoint &b) {return a.FPR<b.FPR;}
};


class Application {
  SequenceSet SS; 
  int trackID, numFiles, windowSize, peakSize, jobID;
  int TP, FP, TN, FN;
  int wmmTP, wmmFP, wmmTN, wmmFN;
  float wmmMin, wmmMax, minLLR, maxLLR;
  bool wantWMM;
  WMM *wmm;
  Array1D< BOOM::Vector<Hit> > windowScores; // [seqID][pos]
  Array1D< BOOM::Vector<Hit> > wmmScores; // [seqID][pos]
  Array2D< BOOM::Vector<Hit> > slices; // [seqID][thresholdID][]
  BOOM::Vector<RocPoint> ROC, wmmROC;
  void getPeaks(EmissionSequence &,BOOM::Vector<Hit> &slices);
  void analyzePeaks(CommandLine &cmd);
  void getSlices(EmissionSequence &,float threshold,BOOM::Vector<Hit> &);
  void resetHits(BOOM::Vector<Hit> &);
  int countOverlaps(BOOM::Vector<Hit> &);
  void getExtrema(BOOM::Vector<Hit> &,float &,float &);
  void sortROC();
  void update(BOOM::Vector<Hit> &windowScores,
	      BOOM::Vector<Hit> &slices,float llrThreshold,
	      int &TP,int &FP,int &TN,int &FN);
  double getAUC(BOOM::Vector<RocPoint> &rocCurve);
  RocPoint *maxPoint(BOOM::Vector<RocPoint> &);
  void getSites(Array1D< BOOM::Vector<Hit> > &,float threshold,
		String outfile);
  void localOptimality(BOOM::Vector<Hit> &windows,float threshold);
public:
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



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"w:b:r:p");
  if(cmd.numArgs()!=10)
    throw string(
"\n\
roc-chain <*.hmm> <fg-state> <bg-state> <site-size> <peak-size> <fastb-dir> <track-name> <min-signal> <max-signal> <jobID>\n\
    where:  -w file = also construct roc curve for PWM in file\n\
            -b file = use background chain in file for PWM\n\
            -r file = write curve into file\n\
            -p = evaluate via peaks rather than uniform slices\n\
"
);
  String hmmFile=cmd.arg(0);
  int fgState=cmd.arg(1).asInt();
  int bgState=cmd.arg(2).asInt();
  windowSize=cmd.arg(3).asInt();
  peakSize=cmd.arg(4).asInt();
  String dir=cmd.arg(5);
  String trackName=cmd.arg(6);
  float minSignal=cmd.arg(7).asFloat();
  float maxSignal=cmd.arg(8).asFloat();
  jobID=cmd.arg(9).asInt();
  wantWMM=cmd.option('w');
  bool wantPeaks=cmd.option('p');

  // Load WMM, if necessary
  GarbageIgnorer gc;
  if(wantWMM) {
    String filename=cmd.optParm('w');
    wmm=new WMM(gc,filename);
  }
  else wmm=NULL;

  // Load hmm
  HMM hmm(hmmFile);
  Schema schema=hmm.getSchema();
  hmm.dropContinuousTracks();
  int numStates=hmm.countStates();
  EmissionLoader loader(schema);
  trackID=schema.lookupContinuousID(trackName);
  const int order=hmm.getOrder();

  // Load sequences
  SS.load(dir,schema);
  numFiles=SS.size();
  windowScores.resize(numFiles);
  slices.resize(numFiles,wantPeaks ? 1 : SIGNAL_SLICES);
  if(wantWMM) wmmScores.resize(numFiles);
  float sliceDelta=(maxSignal-minSignal)/(SIGNAL_SLICES+1);

  // Preprocess sequences
  wmmMin=POSITIVE_INFINITY; wmmMax=NEGATIVE_INFINITY;
  minLLR=POSITIVE_INFINITY; maxLLR=NEGATIVE_INFINITY;
  for(int k=0 ; k<numFiles ; ++k) {
    EmissionSequence &seq=*SS[k];
    const int seqLen=seq.length();
    BOOM::Vector<Hit> &windows=windowScores[k];
    int lastPos=seqLen-windowSize;
    for(int pos=0 ; pos<=lastPos ; ++pos) {
      EmissionSequence *subseq=seq.getSubsequence(pos,windowSize);
      double wmmScore;
      if(wantWMM) {
	Sequence *S=subseq->getDiscreteSeq(0);
	String str;
	S->toString(alphabet,0,S->getLength(),str);
	wmmScore=wmm->getLogP(*S,str,0);
	/*
	Sequence *R=S->reverseComplement(alphabet);
	R->toString(alphabet,0,R->getLength(),str);
	float rScore=wmm->getLogP(*R,str,0);
	if(rScore>wmmScore) wmmScore=rScore;
	delete R;
	*/
	delete S; 
      }
      subseq->recode(order+1);
      double fgScore=hmm.getEmissionProb(fgState,*subseq,true);
      double bgScore=hmm.getEmissionProb(bgState,*subseq,true);
      float score=fgScore-bgScore;
      windows.push_back(Hit(score,pos,pos+windowSize));
      if(score<minLLR) minLLR=score;
      if(score>maxLLR) maxLLR=score;
      if(wantWMM) {
	double wmmLLR=wmmScore-bgScore;
	//double wmmLLR=wmmScore;
	wmmScores[k].push_back(Hit(wmmLLR,pos,pos+windowSize));
	if(wmmLLR<wmmMin) wmmMin=wmmLLR;
	if(wmmLLR>wmmMax) wmmMax=wmmLLR;
      }
      delete subseq;
    }
    if(wantPeaks)
      getPeaks(seq,slices[k][0]);
    else for(int i=0 ; i<SIGNAL_SLICES ; ++i) {
	float threshold=minSignal+(i+1)*sliceDelta;
	getSlices(seq,threshold,slices[k][i]);
      }
  }
  //wmmMin=minLLR; wmmMax=maxLLR;
  //const SMALLEST_MIN=-20;
  //if(minLLR<SMALLEST_MIN) minLLR=SMALLEST_MIN; // ###
  //if(wmmMin<SMALLEST_MIN) wmmMin=SMALLEST_MIN; // ###
  //cerr<<"wmm min="<<wmmMin<<", wmm max="<<wmmMax<<endl;
  //cerr<<"llr min="<<minLLR<<", llr max="<<maxLLR<<endl;

  if(wantPeaks) analyzePeaks(cmd);
  else {
    // Iterate through slicing thresholds
    RocPoint zeroZero, oneOne;
    zeroZero.Sn=zeroZero.FPR=0;
    oneOne.Sn=oneOne.FPR=1;
    for(int i=0 ; i<SIGNAL_SLICES ; ++i) {
      const float sliceThreshold=minSignal+(i+1)*sliceDelta;
      
      // Iterate through LLR thresholds
      ROC.clear(); wmmROC.clear();
      ROC.push_back(oneOne); wmmROC.push_back(oneOne);
      float llrDelta=(maxLLR-minLLR)/(LLR_SLICES+1);
      float wmmDelta=(wmmMax-wmmMin)/(LLR_SLICES+1);
      for(int j=0 ; j<LLR_SLICES ; ++j) {
	const float llrThreshold=minLLR+(j+1)*llrDelta;
	const float wmmThreshold=wmmMin+(j+1)*wmmDelta;
	TP=FP=TN=FN=wmmTP=wmmFP=wmmTN=wmmFN=0;
	
	// Iterate through files
	for(int k=0 ; k<numFiles ; ++k) {
	  EmissionSequence &seq=*SS[k];
	  const int seqLen=seq.length();
	  BOOM::Vector<Hit> &slices=this->slices[k][i];
	  BOOM::Vector<Hit> &windowScores=this->windowScores[k];
	  update(windowScores,slices,llrThreshold,TP,FP,TN,FN);
	  update(wmmScores[k],slices,wmmThreshold,wmmTP,wmmFP,wmmTN,wmmFN);
	}
	RocPoint pt;
	pt.init(TP,FP,FN,TN,llrThreshold);
	ROC.push_back(pt);
	if(wantWMM) { 
	  pt.init(wmmTP,wmmFP,wmmFN,wmmTN,wmmThreshold);
	  wmmROC.push_back(pt);
	}
      }
      ROC.push_back(zeroZero);
      double auc=getAUC(ROC);
      if(wantWMM) {
	wmmROC.push_back(zeroZero);
	double wmmAUC=getAUC(wmmROC);
	cout<<auc<<"\t"<<wmmAUC<<endl;
	if(cmd.option('r')) {
	  const String filename=String("wmm")+i;
	  ofstream os(filename.c_str());
	  int n=wmmROC.size();
	  for(int i=0 ; i<n ; ++i) {
	    //if(i>0 && i%LLR_SLICES==0) os<<endl;
	    RocPoint &pt=wmmROC[i];
	    if(i>0) {
	      RocPoint &r=wmmROC[i-1];
	      if(pt.FPR==r.FPR && pt.Sn==r.Sn) continue;
	    }
	    os<<pt.FPR<<"\t"<<pt.Sn<<endl;
	  }
	}
      }
      else cout<<auc<<endl;
      if(cmd.option('r')) {
	const String filename=cmd.optParm('r')+i;
	ofstream os(filename.c_str());
	int n=ROC.size();
	for(int i=0 ; i<n ; ++i) {
	  //if(i>0 && i%LLR_SLICES==0) os<<endl;
	  RocPoint &pt=ROC[i];
	  if(i>0) {
	    RocPoint &r=ROC[i-1];
	    if(pt.FPR==r.FPR && pt.Sn==r.Sn) continue;
	  }
	  os<<pt.FPR<<"\t"<<pt.Sn<<endl;
	}
      }
    }
  }

  if(wantPeaks && wantWMM) {
    RocPoint *maxChain=maxPoint(ROC), *maxWMM=maxPoint(wmmROC);
    float chainThreshold=maxChain->threshold;
    float wmmThreshold=maxWMM->threshold;

    //wmmThreshold=-20;//### DEBUGGING

    //cout<<"thresholds: "<<chainThreshold<<" "<<wmmThreshold<<endl;
    //cout<<"F: "<<maxChain->F<<" "<<maxWMM->F<<endl;
    getSites(windowScores,chainThreshold,
	     String("chain-sites-")+jobID+".fasta");
    getSites(wmmScores,wmmThreshold,
	     String("pwm-sites-")+jobID+".fasta");
  }

  return 0;
}



void Application::getSites(Array1D< BOOM::Vector<Hit> > &allHits,
			   float threshold,String outfile)
{
  ofstream os(outfile.c_str());
  int numSeqs=SS.size();
  int nextID=1;
  for(int i=0 ; i<numSeqs ; ++i) {
    EmissionSequence &seq=*SS[i];
    BOOM::Vector<Hit> &hits=allHits[i];
    int L=hits.size();
    for(int pos=0 ; pos<L ; ++pos) {
      Hit &hit=hits[pos];
      if(hit.score>=threshold) {
	os<<">"<<nextID++<<endl;
	EmissionSequence *subseq=seq.getSubsequence(hit.begin,
						    hit.end-hit.begin);
	Sequence *S=subseq->getDiscreteSeq(0);
	String &str=*S->toString(alphabet);
	os<<str<<endl;
	delete S; delete subseq; delete &str;
      }
    }
  }
}



RocPoint *Application::maxPoint(BOOM::Vector<RocPoint> &ROC) {
  int n=ROC.size();
  RocPoint *best=&ROC[0];
  for(int i=1 ; i<n ; ++i)
    if(ROC[i].F>best->F) best=&ROC[i];
  return best;
}



void Application::localOptimality(BOOM::Vector<Hit> &windows,float threshold)
{
  int begin;
  bool above=false;
  const int N=windows.size();
  for(int i=0 ; i<N ; ++i) {
    Hit &hit=windows[i];
    if(hit.score>=threshold) {
      if(!above) {
	above=true;
	begin=i;
      }
    }
    else if(above) {
      above=false;
      float max=NEGATIVE_INFINITY;
      int maxPos;
      for(int j=begin ; j<i ; ++j) {
	float val=windows[j].score;
	if(val>max) { max=val; maxPos=j; }
      }
      windows[maxPos].optimal=true;
    }
  }
}



void Application::update(BOOM::Vector<Hit> &wnd,
                         BOOM::Vector<Hit> &slices,float llrThreshold,
                         int &TP,int &FP,int &TN,int &FN)
{
  BOOM::Vector<Hit> windows=wnd;
  int numSlices=slices.size(), numWindows=windows.size();
  resetHits(windows);
  resetHits(slices);
  //localOptimality(windows,llrThreshold);
  int predictions=0;
  for(int k=0 ; k<numWindows ; ++k) 
    if(windows[k].score>=llrThreshold) ++predictions;
  //if(windows[k].score>=llrThreshold && windows[k].optimal) ++predictions;
  for(int i=0 ; i<numSlices ; ++i) {
    Hit &peak=slices[i];
    for(int j=0 ; j<numWindows; ++j) {
      Hit &window=windows[j];
      if(window.begin<peak.end && peak.begin<window.end) {
	if(window.score>=llrThreshold) window.hit=peak.hit=true;
	//if(window.score>=llrThreshold && window.optimal) window.hit=peak.hit=true;
	else window.miss=true;
      }
    }
  }
  int misses=0;
  for(int k=0 ; k<numWindows ; ++k) if(windows[k].miss) ++misses;
  int w=countOverlaps(windows);
  int s=countOverlaps(slices);

  // Sensitivity is evaluated in terms of the peak slices:
  TP+=s;
  FN+=numSlices-s;

  // FPR is evaluated in terms of the actual predictions:
  FP+=predictions-w;
  TN+=numWindows-predictions-misses;
}



void Application::getSlices(EmissionSequence &seq,float threshold,
			    BOOM::Vector<Hit> &slices)
{
  int L=seq.length();
  double prevY=NEGATIVE_INFINITY;
  int begin;
  for(int pos=0 ; pos<L ; ++pos) {
    double y=seq[pos][trackID];
    if(y>=threshold && prevY<threshold) begin=pos;
    else if(y<threshold && prevY>=threshold) {
      slices.push_back(Hit(prevY,begin,pos));
    }
    prevY=y;
  }
}



void Application::resetHits(BOOM::Vector<Hit> &hits)
{
  BOOM::Vector<Hit>::iterator cur=hits.begin(), end=hits.end();
  for(; cur!=end ; ++cur) 
    (*cur).reset();
}



int Application::countOverlaps(BOOM::Vector<Hit> &V)
{
  int c=0;
  BOOM::Vector<Hit>::iterator cur=V.begin(), end=V.end();
  for(; cur!=end ; ++cur)
    if((*cur).hit) ++c;
  return c;
}



void Application::getExtrema(BOOM::Vector<Hit> &V,float &minVal,float &maxVal)
{
  minVal=POSITIVE_INFINITY;
  maxVal=NEGATIVE_INFINITY;
  BOOM::Vector<Hit>::iterator cur=V.begin(), end=V.end();
  for(; cur!=end ; ++cur) {
    float s=(*cur).score;
    if(s<minVal) minVal=s;
    if(s>maxVal) maxVal=s;
  }
}



void Application::sortROC()
{
  RocComparator cmp;
  VectorSorter<RocPoint> sorter(ROC,cmp);
  sorter.sortAscendInPlace();
}



double Application::getAUC(BOOM::Vector<RocPoint> &ROC)
{
  double area=0;
  int n=ROC.size();
  for(int i=1 ; i<n ; ++i) {
    RocPoint &pt=ROC[i], prev=ROC[i-1];
    double x1=pt.FPR, y1=pt.Sn, x2=prev.FPR, y2=prev.Sn;
    //double inc=(x2-x1)*(y1+0.5*abs(y2-y1));
    double inc=(x2-x1)*(y1+0.5*(y2-y1));
    area+=inc;
  }
  return area;
}



void Application::analyzePeaks(CommandLine &cmd)
{
  RocPoint zeroZero, oneOne;
  zeroZero.Sn=zeroZero.FPR=0;
  oneOne.Sn=oneOne.FPR=1;
    
  // Iterate through LLR thresholds
  ROC.clear(); wmmROC.clear();
  ROC.push_back(oneOne); wmmROC.push_back(oneOne);
  float llrDelta=(maxLLR-minLLR)/(LLR_SLICES+1);
  float wmmDelta=(wmmMax-wmmMin)/(LLR_SLICES+1);
  for(int j=0 ; j<LLR_SLICES ; ++j) {
    const float llrThreshold=minLLR+(j+1)*llrDelta;
    const float wmmThreshold=wmmMin+(j+1)*wmmDelta;
    TP=FP=TN=FN=wmmTP=wmmFP=wmmTN=wmmFN=0;
    
    // Iterate through files
    for(int k=0 ; k<numFiles ; ++k) {
      EmissionSequence &seq=*SS[k];
      const int seqLen=seq.length();
      BOOM::Vector<Hit> &slices=this->slices[k][0];
      update(windowScores[k],slices,llrThreshold,TP,FP,TN,FN);
      update(wmmScores[k],slices,wmmThreshold,wmmTP,wmmFP,wmmTN,wmmFN);
    }
    RocPoint pt;
    pt.init(TP,FP,FN,TN,llrThreshold);
    ROC.push_back(pt);
    if(wantWMM) { 
      pt.init(wmmTP,wmmFP,wmmFN,wmmTN,wmmThreshold);
      wmmROC.push_back(pt);
    }
  }
  ROC.push_back(zeroZero);
  double auc=getAUC(ROC);
  if(wantWMM) {
    wmmROC.push_back(zeroZero);
    double wmmAUC=getAUC(wmmROC);
    cout<<auc<<"\t"<<wmmAUC<<endl;
    if(cmd.option('r')) {
      const String filename=String("wmm.roc-")+jobID;
      ofstream os(filename.c_str());
      int n=wmmROC.size();
      for(int i=0 ; i<n ; ++i) {
	//if(i>0 && i%LLR_SLICES==0) os<<endl;
	RocPoint &pt=wmmROC[i];
	if(i>0) {
	  RocPoint &r=wmmROC[i-1];
	  if(pt.FPR==r.FPR && pt.Sn==r.Sn) continue;
	}
	os<<pt.FPR<<"\t"<<pt.Sn<<endl;
      }
    }
  }
  else cout<<auc<<endl;
  if(cmd.option('r')) {
    const String filename=cmd.optParm('r')+".roc-"+jobID;
    ofstream os(filename.c_str());
    int n=ROC.size();
    for(int i=0 ; i<n ; ++i) {
      //if(i>0 && i%LLR_SLICES==0) os<<endl;
      RocPoint &pt=ROC[i];
      if(i>0) {
	RocPoint &r=ROC[i-1];
	if(pt.FPR==r.FPR && pt.Sn==r.Sn) continue;
      }
      os<<pt.FPR<<"\t"<<pt.Sn<<endl;
    }
  }
}



void Application::getPeaks(EmissionSequence &seq,
			   BOOM::Vector<Hit> &slices)
{
  int L=seq.length();
  int lastPos=L-peakSize;
  for(int pos=0 ; pos<=lastPos ; ++pos) {
    const int centerPos=pos+peakSize/2;
    const double centerVal=seq[centerPos][trackID];
    const int end=pos+peakSize;
    bool isMax=true;
    for(int i=pos ; i<end ; ++i)
      if(i!=centerPos && seq[i][trackID]>=centerVal) {isMax=false;break;}
    if(isMax) {
      slices.push_back(Hit(centerVal,pos,pos+peakSize));
      //pos+=peakSize;
    }
  }
}
