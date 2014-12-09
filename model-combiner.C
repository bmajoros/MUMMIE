/**************************************************************
 model-combiner.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/

#include <iostream>
#include <fstream>
#include "BOOM/File.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Map.H"
#include "BOOM/Vector.H"
#include "BOOM/HigherOrderAlphabet.H"
#include "HMM.H"
#include "HMMreader.H"
using namespace std;


class Application {
  HMM *metaModel;
  int numSubmodels;
  BOOM::Vector<int> submodelOrigins; // state ids in composite HMM
  BOOM::Vector<String> submodelNames;
  BOOM::Vector<HMM*> *loadSubmodels(const String &filename);
  HMM *merge(BOOM::Vector<HMM*> &submodels);
  void copyIn(HMM &submodel,HMM &compositeHMM,int startingAt);
  void linkModels(BOOM::Vector<HMM*> &,int fromModelId,
		  int toModelId,HMM &compositeHMM,double metaP);
  void closeTransitivity(int fromState,int toModelId,double fromP,
			 HMM &compositeHMM,BOOM::Vector<HMM*> &);
public:
  Application() {}
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
  catch(string s)
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
  throw "model-combiner <meta-model.hmms> <submodels.txt> <outfile.hmm>";
}


int Application::go(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=3) usage();
  String structureFile=cmd.arg(0);
  String submodelListFile=cmd.arg(1);
  String outfile=cmd.arg(2);

  // Load the MetaModel
  HMMreader reader;
  metaModel=reader.read(structureFile);

  // Load the submodels
  BOOM::Vector<HMM*> &submodels=*loadSubmodels(submodelListFile);
  metaModel->getSchema()=submodels[1]->getSchema();
  numSubmodels=submodels.size()-1;
  submodels[0]=NULL;

  // Merge the submodels into a composite model
  HMM *compositeHMM=merge(submodels);
  compositeHMM->normalizeTransitions();

  // Save the results
  compositeHMM->save(outfile);

  return 0;
}



BOOM::Vector<HMM*> *Application::loadSubmodels(const String &filename)
{
  submodelNames.push_back("");
  BOOM::Vector<HMM*> *models=new BOOM::Vector<HMM*>;
  File file(filename);
  while(!file.eof()) {
    String line=file.readLine();
    if(file.eof()) break;
    BOOM::Vector<String> &fields=*line.getFields("= \t\n+,");
    int numFields=fields.size();
    if(numFields==1)
      throw String("Syntax error in submodel list file: ")+line;
    if(numFields>0) {
      String firstField=fields[0];
      int modelId=firstField.asInt();
      String modelFilename=fields[1];
      delete &fields;
      submodelNames.push_back(modelFilename);
      HMM *submodel=new HMM(modelFilename);
      while(modelId>=models->size()) models->push_back(NULL);
      (*models)[modelId]=submodel;
    }
  }
  return models;
}



void Application::copyIn(HMM &submodel,HMM &compositeHMM,int delta)
{
  // Copy transitions
  int numStates=submodel.countStates();
  for(int i=1 ; i<numStates ; ++i) {
    for(int j=1 ; j<numStates ; ++j)
      compositeHMM.setTransitionProb(i+delta,j+delta,
				     submodel.getTransitionProb(i,j));
    compositeHMM.setTransitionProb(i+delta,0,
				   submodel.getTransitionProb(i,0));
  }
  
  // Copy emissions
  for(int i=1 ; i<numStates ; ++i)
    compositeHMM.copyEmissionDistr(submodel,i,i+delta);
}



HMM *Application::merge(BOOM::Vector<HMM*> &submodels)
{
  // Determine how many states we need in the composite model
  int totalStates=1; // has a silent zero state
  int maxOrder=0;
  for(int i=1 ; i<=numSubmodels ; ++i) {
    HMM &submodel=*submodels[i];
    totalStates+=submodel.countStates()-1;//don't count zero states
    int order=submodel.getOrder();
    if(order>maxOrder) maxOrder=order;
  }

  // Allocate a composite model
  HMM &compositeHMM=*new HMM(totalStates);
  Array2D< Array1D<double> > &emitDistr=compositeHMM.getDiscreteEmitDistr();
  HMM *firstSubmodel=submodels[1];
  Schema schema=firstSubmodel->getSchema();
  int numDiscrete=schema.getNumDiscrete();
  emitDistr.resize(totalStates,numDiscrete);
  /*
  for(int i=0 ; i<numDiscrete ; ++i) {
    Alphabet &alpha=schema.getAlphabet(i);
    HigherOrderAlphabet H(alpha,order+1);
    int numNmers=H.getNumNmers();
    for(int q=1 ; q<totalStates ; ++q) emitDistr[q][i].resize(numNmers);//###
  }
  */
  compositeHMM.getSchema()=schema;
  compositeHMM.changeOrder(maxOrder);
  
  // Copy the submodels into distinct, unconnected regions of the
  // composite model
  int nextFreeState=1;
  submodelOrigins.push_back(0);
  for(int modelId=1 ; modelId<=numSubmodels ; ++modelId) {
    HMM &submodel=*submodels[modelId];
    copyIn(submodel,compositeHMM,nextFreeState-1);
    submodelOrigins.push_back(nextFreeState);
    cout << "submodel #" << modelId << " has " 
	 << submodel.countStates()-1 
	 << " states, numbered " << nextFreeState << "-"
	 << nextFreeState+submodel.countStates()-2 
	 << "\t" << submodelNames[modelId] << endl;
    nextFreeState+=submodel.countStates()-1;
  }
  
  // Establish connections between the submodels, as dictated by
  // the MetaModel
  for(int from=1 ; from<=numSubmodels ; ++from)
    for(int to=1 ; to<=numSubmodels ; ++to)
      if(metaModel->doesTransitionExist(from,to))
	linkModels(submodels,from,to,compositeHMM,
		   metaModel->getTransitionProb(from,to));

  // Establish connections from state 0 in the metaModel
  for(int toModelId=1 ; toModelId<=numSubmodels ; ++toModelId)
    if(metaModel->doesTransitionExist(0,toModelId)) {
      double metaP=metaModel->getTransitionProb(0,toModelId);
      HMM &toModel=*submodels[toModelId];
      int toModelDelta=submodelOrigins[toModelId]-1;
      int numToModelStates=toModel.countStates();
      for(int toState=1 ; toState<numToModelStates ; ++toState)	{
	if(toModel.doesTransitionExist(0,toState)) {
	  double newP=metaP*toModel.getTransitionProb(0,toState);
	  compositeHMM.setTransitionProb(0,toState+toModelDelta,newP);
	}
      }

      // Handle 0->0 transitions
      if(toModel.doesTransitionExist(0,0)) {
	double newP=metaP*toModel.getTransitionProb(0,0);
	closeTransitivity(0,toModelId,newP,compositeHMM,submodels);
      }
    }

  // Establish connections to state 0 in the metaModel
  for(int fromModelId=1 ; fromModelId<=numSubmodels ; ++fromModelId)
    if(metaModel->doesTransitionExist(fromModelId,0)) {
      double metaP=metaModel->getTransitionProb(fromModelId,0);
      HMM &fromModel=*submodels[fromModelId];
      int fromModelDelta=submodelOrigins[fromModelId]-1;
      int numFromModelStates=fromModel.countStates();
      for(int fromState=1 ; fromState<numFromModelStates ; ++fromState)	{
	    if(fromModel.doesTransitionExist(fromState,0)) {
	      double newP=metaP*fromModel.getTransitionProb(fromState,0);
	      compositeHMM.setTransitionProb(fromState+fromModelDelta,0,
					     newP);
	    }
      }
    }
    else {
      HMM &fromModel=*submodels[fromModelId];
      int fromModelDelta=submodelOrigins[fromModelId]-1;
      int numFromModelStates=fromModel.countStates();
      for(int fromState=1 ; fromState<numFromModelStates ; ++fromState) {
	compositeHMM.setTransitionProb(fromState+fromModelDelta,0,0);
      }
    }

  return &compositeHMM;
}



void Application::linkModels(BOOM::Vector<HMM*> &submodels,
			     int fromModelId,int toModelId,
			     HMM &compositeHMM,double metaP)
{
  /*
    This method takes every state in the "from" model having a
    transition to zero and redirects that transition so that
    it points to the state(s) that the zero in the "to" model
    can transition to.
   */

  HMM &fromModel=*submodels[fromModelId];
  HMM &toModel=*submodels[toModelId];
  int fromModelDelta=submodelOrigins[fromModelId]-1;
  int toModelDelta=submodelOrigins[toModelId]-1;
  int numToModelStates=toModel.countStates();
  int numFromModelStates=fromModel.countStates();

  for(int fromState=1 ; fromState<numFromModelStates ; ++fromState) {
    if(!fromModel.doesTransitionExist(fromState,0)) continue;
    double fromP=fromModel.getTransitionProb(fromState,0);
    for(int toState=1 ; toState<numToModelStates ; ++toState) {
      if(!toModel.doesTransitionExist(0,toState)) continue;
      double toP=toModel.getTransitionProb(0,toState);
      double newP=fromP*metaP*toP;
      if(compositeHMM.doesTransitionExist(fromState+fromModelDelta,
					  toState+toModelDelta)) {
	cout<<"fromModelDelta="<<fromModelDelta<<" toModelDelta="<<toModelDelta<<" numFromModelStates="<<numFromModelStates<<" numToModelStates="<<numToModelStates<<" fromModelId="<<fromModelId<<" toModelId="<<toModelId<<" fromState="<<fromState<<" toState="<<toState<<" transP="<<compositeHMM.getTransitionProb(fromState+fromModelDelta,toState+toModelDelta)<<" newP="<<newP<<endl;
	throw String("ERROR: transition already exists!");
      }
      compositeHMM.setTransitionProb(fromState+fromModelDelta,
				     toState+toModelDelta,
				     newP);
    }
    if(toModel.doesTransitionExist(0,0)) {
      double newP=fromP*metaP*toModel.getTransitionProb(0,0);
      closeTransitivity(fromState+fromModelDelta,
			toModelId,newP,compositeHMM,submodels);
    }
  }
}



void Application::closeTransitivity(int fromState,int intermediateModelId,
				    double fromP,HMM &compositeHMM,
				    BOOM::Vector<HMM*> &submodels)
{
  for(int toModelId=1 ; toModelId<=numSubmodels ; ++toModelId)
    if(metaModel->doesTransitionExist(intermediateModelId,toModelId)) {
      double metaP=
	metaModel->getTransitionProb(intermediateModelId,toModelId);
      HMM &toModel=*submodels[toModelId];
      int toModelDelta=submodelOrigins[toModelId]-1;
      int numToModelStates=toModel.countStates();
      for(int toState=1 ; toState<numToModelStates ; ++toState) {
	if(!toModel.doesTransitionExist(0,toState)) continue;
	double toP=toModel.getTransitionProb(0,toState);
	double newP=fromP*metaP*toP;
	if(compositeHMM.doesTransitionExist(fromState,
					    toState+toModelDelta))
	  throw String("Error: transition already exists!");
	compositeHMM.setTransitionProb(fromState,
				       toState+toModelDelta,
				       newP);
      }
      if(toModel.doesTransitionExist(0,0)) {
	double newP=fromP*metaP*toModel.getTransitionProb(0,0);
	closeTransitivity(fromState,toModelId,newP,compositeHMM,
			  submodels);
      }
    }
}




