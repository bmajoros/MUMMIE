/*
 roc.C
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
#include "BOOM/VectorSorter.H"
#include "BOOM/PureDnaAlphabet.H"
#include "BOOM/Map.H"
#include "BOOM/GffReader.H"
#include "EGGS/WMM.H"
#include "EGGS/GarbageCollector.H"
#include "FastViterbi.H"
#include "EmissionLoader.H"
#include "SequenceSet.H"
using namespace std;

BOOM::Alphabet alphabet=PureDnaAlphabet::global();

const int SIGNAL_SLICES=100; //50;

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
  int trackID, numFiles;
  int TP, FP, TN, FN;
  Map<String,int> substrateToIndex;
  Regex substrateRegex;
  Array1D< BOOM::Vector<Hit> > predictions; // [seqID][]
  Array2D< BOOM::Vector<Hit> > slices; // [seqID][thresholdID][]
  BOOM::Vector<RocPoint> ROC;
  void getSlices(EmissionSequence &,float threshold,BOOM::Vector<Hit> &);
  void resetHits(BOOM::Vector<Hit> &);
  int countOverlaps(BOOM::Vector<Hit> &);
  void getExtrema(BOOM::Vector<Hit> &,float &,float &);
  void sortROC();
  void update(BOOM::Vector<Hit> &predictions,
	      BOOM::Vector<Hit> &slices,
	      int &TP,int &FP,int &TN,int &FN);
  double getAUC(BOOM::Vector<RocPoint> &rocCurve);
  RocPoint *maxPoint(BOOM::Vector<RocPoint> &);
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
  : substrateRegex("([^/]+)\\.fastb")
{
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=7)
    throw string(
"\n\
roc <in.gff> <out.roc> <fastb-dir> <track-name> <min-signal> <max-signal> <schema>\n\
"
);
  String gffFile=cmd.arg(0);
  String rocFile=cmd.arg(1);
  String dir=cmd.arg(2);
  String trackName=cmd.arg(3);
  float minSignal=cmd.arg(4).asFloat();
  float maxSignal=cmd.arg(5).asFloat();
  String schemaFile=cmd.arg(6);

  // Load sequences
  Schema schema(schemaFile);
  trackID=schema.lookupContinuousID(trackName);
  SS.load(dir,schema);
  numFiles=SS.size();
  predictions.resize(numFiles);
  slices.resize(numFiles,SIGNAL_SLICES);
  float sliceDelta=(maxSignal-minSignal)/(SIGNAL_SLICES+1);
  for(int i=0 ; i<numFiles ; ++i) {
    EmissionSequence &seq=*SS[i];
    const String &filename=seq.getFilename();
    if(!substrateRegex.search(filename))
      throw String("can't parse filename: ")+filename;
    const String substrate=substrateRegex[1];
    substrateToIndex[substrate]=i;
  }

  // Load predictions
  GffReader reader(gffFile);
  while(1) {
    GffFeature *feat=reader.nextFeature();
    if(!feat) break;
    const String &substrate=feat->getSubstrate();
    if(!substrateToIndex.isDefined(substrate)) {
      //throw String("undefined substrate: ")+substrate;
      delete feat;
      continue;
    }
    int index=substrateToIndex[substrate];
    int begin=feat->getBegin(), end=feat->getEnd();
    float score=feat->getScore();
    predictions[index].push_back(Hit(score,begin,end));
    delete feat;
  }

  // Preprocess sequences
  for(int k=0 ; k<numFiles ; ++k) {
    EmissionSequence &seq=*SS[k];
    const int seqLen=seq.length();
    for(int i=0 ; i<SIGNAL_SLICES ; ++i) {
      float threshold=minSignal+(i+1)*sliceDelta;
      getSlices(seq,threshold,slices[k][i]);
    }
  }

  // Iterate through slicing thresholds
  RocPoint zeroZero, oneOne;
  zeroZero.Sn=zeroZero.FPR=0;
  oneOne.Sn=oneOne.FPR=1;
  for(int i=0 ; i<SIGNAL_SLICES ; ++i) {
    const float sliceThreshold=minSignal+(i+1)*sliceDelta;

    ROC.clear();
    ROC.push_back(oneOne);
    TP=FP=TN=FN=0;
    
    // Iterate through files
    for(int k=0 ; k<numFiles ; ++k) {
      EmissionSequence &seq=*SS[k];
      const int seqLen=seq.length();
      BOOM::Vector<Hit> &slices=this->slices[k][i];
      BOOM::Vector<Hit> &windowScores=this->predictions[k];
      update(windowScores,slices,TP,FP,TN,FN);
    }
    RocPoint pt;
    pt.init(TP,FP,FN,TN,0);
    ROC.push_back(pt);
  }
  ROC.push_back(zeroZero);
  double auc=getAUC(ROC);
  cout<<auc<<endl;
  
  ofstream os(rocFile.c_str());
  int n=ROC.size();
  for(int i=0 ; i<n ; ++i) {
    RocPoint &pt=ROC[i];
    if(i>0) {
      RocPoint &r=ROC[i-1];
      if(pt.FPR==r.FPR && pt.Sn==r.Sn) continue;
    }
    os<<pt.FPR<<"\t"<<pt.Sn<<endl;
  }

  return 0;
}



RocPoint *Application::maxPoint(BOOM::Vector<RocPoint> &ROC) {
  int n=ROC.size();
  RocPoint *best=&ROC[0];
  for(int i=1 ; i<n ; ++i)
    if(ROC[i].F>best->F) best=&ROC[i];
  return best;
}



void Application::update(BOOM::Vector<Hit> &predictionsV,
			 BOOM::Vector<Hit> &slices,
			 int &TP,int &FP,int &TN,int &FN)
{
  BOOM::Vector<Hit> windows=predictionsV;
  int numSlices=slices.size(), numWindows=windows.size();
  resetHits(windows);
  resetHits(slices);
  int predictions=windows.size();
  for(int i=0 ; i<numSlices ; ++i) {
    Hit &peak=slices[i];
    for(int j=0 ; j<predictions; ++j) {
      Hit &window=windows[j];
      if(window.begin<peak.end && peak.begin<window.end)
	window.hit=peak.hit=true;
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
    double inc=(x2-x1)*(y1+0.5*(y2-y1));
    area+=inc;
  }
  return area;
}
