/****************************************************************
 test-multigauss.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Array1D.H"
#include "BOOM/Array3D.H"
#include "BOOM/GSL/MultiGauss.H"
#include "BOOM/GSL/GaussianDistribution.H"
using namespace std;
using namespace BOOM;

const double MIN=2, MAX=8;
const int NUM_CELLS=10;
const int NUM_POINTS=100000;

class Application
{
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
  {
    // ctor
  }



int discretize(double d)
{
  const double range=MAX-MIN;
  const double cellSize=range/NUM_CELLS;
  int x=(d-MIN)/cellSize;
  //if(x<0) cout<<x<<" XXX"<<endl;
  if(x<0) x=0;
  //if(x>=NUM_CELLS) cout<<x<<" YYY"<<endl;
  if(x>=NUM_CELLS) x=NUM_CELLS-1;
  return x;
}


int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=0)
    throw String("test-multigauss <no arguments>");

  GSL::Vector means(3);
  means[0]=(MIN+MAX)/2.0;//0.5;
  means[1]=(MIN+MAX)/2.0;//5.0;
  means[2]=(MIN+MAX)/2.0;//100;

  GSL::Matrix cov(3,3);
  cov(0,0)=cov(1,1)=cov(2,2)=1;
  //cov(2,2)=50;
  cov(0,1)=cov(1,0)=0;
  cov(0,2)=cov(2,0)=0;
  cov(1,2)=cov(2,1)=0;

  GSL::MultiGauss G(means,cov);

  //GSL::GaussianDistribution g(3,1);
  //ofstream os1("test1.txt"), os2("test2.txt"), os3("test3.txt");
  //ofstream os12("test12.txt"), os23("test23.txt"), os13("test13.txt");
  //Array1D<double> E(3);
  //E.setAllTo(0);
  GSL::Vector x(3);

  /*
  for(x[0]=2 ; x[0]<4 ; x[0]+=0.01) {
    //E[0]+=g.probabilityEquals(x[0])*x[0]*0.01;
    for(x[1]=2 ; x[1]<4 ; x[1]+=0.01) {
      for(x[2]=2 ; x[2]<4 ; x[2]+=0.01) {
	double pdf=G.density(x);
	//cout<<(fabs(x[2]-x[1])+fabs(x[2]-x[0])+fabs(x[1]-x[0]))/3<<"\t"<<pdf<<endl;
	for(int j=0 ; j<3 ; ++j)
	  
	  //E[j]+=pdf*x[j];
	  //cout<<x<<"\t"<<g.probabilityEquals(x[j])<<endl;
	  ;//E[j]+=g.probabilityEquals(x[j])*x[j]/2e6;
      }
    }
  }
  cout<<E<<endl;
  return 0;
  */

  Array3D<double> counts(NUM_CELLS,NUM_CELLS,NUM_CELLS);
  counts.setAllTo(0.0);
  for(int i=0 ; i<NUM_POINTS ; ++i) {
    GSL::Vector x;
    G.sample(x);
    /*
    os1<<x[0]<<endl;
    os2<<x[1]<<endl;
    os3<<x[2]<<endl;
    os12<<x[0]<<"\t"<<x[1]<<endl;
    os23<<x[1]<<"\t"<<x[2]<<endl;
    os13<<x[0]<<"\t"<<x[2]<<endl;
    */
    //double pdf=G.density(x);
    int y[3];
    for(int j=0;j<3;++j) y[j]=discretize(x[j]);
    counts[y[0]][y[1]][y[2]]+=1.0;///NUM_POINTS;

    //for(int j=0 ; j<3 ; ++j) E[j]+=g.probabilityEquals(x[j])*x[j];
    //cout<<x[2]<<"\t"<<pdf<<endl;
    //cout<<(fabs(x[2]-x[1])+fabs(x[2]-x[0])+fabs(x[1]-x[0]))/3<<"\t"<<pdf<<endl;
    //cout<<pdf<<endl;
    //cout<<g.probabilityEquals(x[0])<<endl;
    //double z=g.random();
    //cout<<z<<"\t"<<g.probabilityEquals(z)<<endl;
    //cout<<g.probabilityEquals(z)<<endl;
  }

  double pdfSum=0, countSum=0;
  double cellSize=(MAX-MIN)/NUM_CELLS;
  for(x[0]=MIN+cellSize/2 ; x[0]<MAX ; x[0]+=cellSize) {
    for(x[1]=MIN+cellSize/2 ; x[1]<MAX ; x[1]+=cellSize) {
      for(x[2]=MIN+cellSize/2 ; x[2]<MAX ; x[2]+=cellSize) {
	double pdf=G.density(x);
	pdfSum+=pdf*cellSize*cellSize*cellSize;
      }
    }
  }
  for(x[0]=MIN+cellSize/2 ; x[0]<MAX ; x[0]+=cellSize) {
    for(x[1]=MIN+cellSize/2 ; x[1]<MAX ; x[1]+=cellSize) {
      for(x[2]=MIN+cellSize/2 ; x[2]<MAX ; x[2]+=cellSize) {
	double pdf=G.density(x);
	if(pdf>1.00001) throw String("pdf>1: ")+pdf;
	double c=counts[discretize(x[0])][discretize(x[1])][discretize(x[2])]
	  /NUM_POINTS;
	cout<<c<<"\t"<<pdf*cellSize*cellSize*cellSize<<endl;
	countSum+=c;
      }
    }
  }
  cerr<<"ZZZ "<<pdfSum<<"\t"<<countSum<<endl;

  //cout<<E<<endl;
  return 0;
}

