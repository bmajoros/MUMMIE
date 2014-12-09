/**************************************************************
 random-HMM.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#include <string>
#include <iostream>
#include "BOOM/Random.H"
#include "BOOM/DnaAlphabet.H"
#include "BOOM/CommandLine.H"
#include "HMMbuilder.H"
#include "Schema.H"
using namespace std;

class Application
{
  DnaAlphabet alphabet;
  int numStates;
  bool wantUniformChains;
  TransitionList *getTemplate(String typeName);
  TransitionList *constructList(int *arrayOfPairs);
  TransitionList *type0();
  TransitionList *typeZ();
  TransitionList *typeI();
  TransitionList *typeII();
  TransitionList *typeIII();
  TransitionList *typeIV();
  TransitionList *typeV();
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
  CommandLine cmd(argc,argv,"t:s:c:u");
  if(cmd.numArgs()!=6)
    throw "\n\
random-HMM [options] <num-states> <connectivity[0-1]> <num-mixture-components> <*.schema> <order> <outfile>\n\
  where: -t *.hmm : use given hmm as template for transition constraints\n\
         -s seed : use given randomization seed\n\
         -c X : use canonical topology of type X (I, II, III, IV, ...)\n\
         -u : use uniform distribution for discrete chains\n\
";
  numStates=cmd.arg(0).asInt();
  double connectivity=cmd.arg(1).asDouble();
  int numComp=cmd.arg(2).asInt();
  String schemaFile=cmd.arg(3);
  int order=cmd.arg(4).asInt();
  String outfile=cmd.arg(5);
  if(cmd.option('s')) SeedRandomizer(cmd.optParm('s').asUnsigned());
  else randomize();
  wantUniformChains=cmd.option('u');

  // Generate random HMM
  Schema schema(schemaFile);
  HMMbuilder builder;
  TransitionList *transList=NULL;
  if(cmd.option('c')) transList=getTemplate(cmd.optParm('c'));
  HMM *hmm=builder.randomHMM(numStates,connectivity,numComp,schema,order,
			     transList,wantUniformChains);
  if(cmd.option('c') && cmd.optParm('c')=="V")
    for(int q=1 ; q<numStates ; ++q)
      hmm->changeStateOrder(q,q-1);

  // Enforce transition constraints, if provided
  if(cmd.option('t')) {
    String filename=cmd.optParm('t');
    HMM constraints(filename);
    for(int i=0 ; i<numStates ; ++i) {
      for(int j=0 ; j<numStates ; ++j) {
	if(!constraints.doesTransitionExist(i,j))
	  hmm->setTransitionProb(i,j,0.0);
      }
    }
    hmm->normalizeTransitions();
  }

  // Save the results
  hmm->save(outfile);
  return 0;
}



TransitionList *Application::getTemplate(String typeName)
{
  if(typeName=="Z") return typeZ();
  if(typeName=="0") return type0();
  if(typeName=="I") return typeI();
  if(typeName=="II") return typeII();
  if(typeName=="III") return typeIII();
  if(typeName=="IV") return typeIV();
  if(typeName=="V") return typeV();
}



TransitionList *Application::constructList(int *A)
{
  TransitionList *L=new TransitionList;
  for(int i=0 ; A[i]>=0 ; i+=2)
    L->push_back(StatePair(A[i],A[i+1]));
  return L;
}



TransitionList *Application::typeI() 
{
  int pairs[]={0,1, 1,0, 1,1, 1,2, 2,1, 2,2, -1};
  return constructList(pairs);
}



TransitionList *Application::typeII()
{
  int pairs[]={0,1, 1,0, 1,1, 1,2, 2,1, 2,2, 2,3, 3,2, 3,3, -1};
  return constructList(pairs);
}



TransitionList *Application::typeIII()
{
  int pairs[]={0,1, 1,0, 1,1, 1,2, 2,2, 2,3, 3,3, 3,4, 4,2, 4,4, 4,1, -1};
  return constructList(pairs);
}



TransitionList *Application::typeIV()
{
  int pairs[]={0,1, 1,0, 1,1, 1,2, 2,2, 2,3, 3,3, 3,4, 4,4, 4,5,
	       5,3, 5,5, 5,6, 6,2, 6,6, 6,1, -1};
  return constructList(pairs);
}



TransitionList *Application::typeV()
{
  TransitionList *t=new TransitionList;
  for(int q=0 ; q<numStates-1 ; ++q) t->push_back(StatePair(q,q+1));
  t->push_back(StatePair(numStates-1,0));
  return t;
}



TransitionList *Application::type0()
{
  int pairs[]={0,1, 1,0, 1,1, -1};
  return constructList(pairs);
}



TransitionList *Application::typeZ()
{
  int pairs[]={0,1, 1,0, -1};
  return constructList(pairs);
}




