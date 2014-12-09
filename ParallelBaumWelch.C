/**************************************************************
 ParallelBaumWelch.C : the Baum-Welch Expectation Maximization algorithm
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/

#include "ParallelBaumWelch.H"

BWThread::BWThread(int states,int numAlpha,HMM &hmm,
                   Vector<Sequence*> &trainingSet,int first,int last)
    : E(states,numAlpha),
      A(states,states),
      hmm(hmm),
      trainingSet(trainingSet),
      first(first),
      last(last),
      likelihoods(trainingSet.size()),
      scalingRatios(1),
      states(states),
      alphabetSize(numAlpha),
      likelihoodCorrection(0)
{
}



void BWThread::f()
{
    A.setAllTo(0); // don't use pseudocounts--they cause likelihood drops!
    E.setAllTo(0); // (unless you re-do the math to account for them...)
    //for(int symbol=0 ; symbol<alphabetSize ; ++symbol) E[0][symbol]=0;
    
    for(int j=first ; j<=last ; ++j)
    {
        Sequence &sequence=*trainingSet[j];
        int sequenceLength=sequence.getLength();
        if(sequenceLength>0)
        {
            /***********************************************************
		Run the Forward algorithm to compute f(k,i)=P[a model M in
		the initial state will emit sequence x(1)..x(i) and reside
		in state k when emitting xi at time i]:
            ***********************************************************/
            ForwardAlgorithm f(hmm,sequence);
            const Array1D<double> &scalingValues=
		f.getScalingFactors();
            likelihoods[j]=f.getScaledP();
            for(int i=1 ; i<=sequenceLength ; ++i)
		likelihoodCorrection+=log(scalingValues[i]);
            
            /***********************************************************
		Run the Backward algorithm to compute b(k,i)=P[M will 
		next emit the string x(i+1)..x(L) and then terminate |
		M is currently in state k]:
            ***********************************************************/
            BackwardAlgorithm b(hmm,sequence);
            
            /***********************************************************
		Compute scaling ratios so the scaling values can be
		properly cancelled during the updating of the A and E
		counts:
            ***********************************************************/
            scalingRatios.resize(sequenceLength+1);
            scalingRatios[sequenceLength]=
		1/f.getScalingFactors()[sequenceLength];
            for(int i=sequenceLength-1 ; i>0 ; --i)
            {
                double r=scalingRatios[i+1] *
		    b.getScalingFactors()[i+1]/f.getScalingFactors()[i];
                scalingRatios[i]=r;
                if(isinf(r) || isnan(r) || r==0)
                {
		    cerr << "WARNING: scaling ratio "<<i<<"="
			 <<scalingRatios[i]<<endl;
                    cout<<scalingRatios[i+1]<<" "<<f.getScalingFactors()[i]<<" "<<b.getScalingFactors()[i+1]<<endl;
                    //cout<<alphabet.lookup(sequence[i-1])<<","<<alphabet.lookup(sequence[i])<<","<<alphabet.lookup(sequence[i+1])<<endl;
                    throw "quitting";
                }
            }
            
            /***********************************************************
		Update the expected number of transitions and emissions
		used in generating this sequence:
            ***********************************************************/
            reviseExpectedTransCounts(sequenceLength,f,b,sequence);
            reviseExpectedEmitCounts(sequenceLength,f,b,sequence);
        }
    }
}



ParallelBaumWelch::ParallelBaumWelch(HMM &hmm,long maxIterations,
                                     Vector<Sequence*> &trainingSet,
                                     int numThreads)
  : hmm(hmm), 
    numHmmStates(hmm.countStates()),alphabet(hmm.getAlphabet()),
    numTrain(trainingSet.size()),likelihoods(trainingSet.size()),
    trainingSet(trainingSet),
    alphabetSize(alphabet.getNumElements()),scalingRatios(1),
    maxIterations(maxIterations),
    E(hmm.countStates(),alphabet.getNumElements()),
    A(hmm.countStates(),hmm.countStates()),
    numThreads(numThreads)
{
  mainAlgorithm();
}



/***************************************************************
  This is the main Baum-Welch Expectation Maximization algorithm
***************************************************************/
void ParallelBaumWelch::mainAlgorithm()
{
  double logLikelihood, deltaLikelihood=logLikelihoodThreshold+1, 
    oldLogLikelihood=0;

  /*************************************************************
    Various initialization steps
  *************************************************************/
  int numEmptyStrings=0;
  for(int j=0 ; j<numTrain ; ++j)
    {
      Sequence &sequence=*trainingSet[j];
      if(sequence.getLength()==0) ++numEmptyStrings;
    }
  bool allowEmptyStrings=hmm.doesTransitionExist(0,0);
  hmm.setTransitionProb(0,0,0);
  progress.start(maxIterations);

  /***********************************************************
    Main loop
  ***********************************************************/
  for(int i=0 ; i<maxIterations ; ++i)
    {
      A.setAllTo(0); // don't use pseudocounts--they cause likelihood drops!
      E.setAllTo(0); // (unless you re-do the math to account for them...)
      for(int symbol=0 ; symbol<alphabetSize ; ++symbol) E[0][symbol]=0;
      
      // Spawn threads
      float numPerThread=numTrain/float(numThreads);
      float begin=0;
      int first, last;
      Array1D<BWThread*> threads(numThreads);
      for(int i=0 ; i<numThreads ; ++i)
      {
          float end=begin+numPerThread-1;
          first=int(begin);
          last=int(end);
          if(i==numThreads-1) last=numTrain-1;
          if(last>=numTrain) last=numTrain-1;
          BWThread *thread=new BWThread(numHmmStates,alphabetSize,hmm,
                                        trainingSet,first,last);
          thread->start();
          begin=end;
          threads[i]=thread;
      }

      // Recombine answers from different threads
      A.setAllTo(0); // don't use pseudocounts--they cause likelihood drops!
      E.setAllTo(0); // (unless you re-do the math to account for them...)
      likelihoods.setAllTo(0);
      double likelihoodCorrection=0;
      for(int i=0 ; i<numThreads ; ++i)
      {
          BWThread *thread=threads[i];
          thread->join();
          for(int i=0 ; i<numHmmStates ; ++i)
              for(int j=0 ; j<numHmmStates ; ++j)
                  A[i][j]+=thread->A[i][j];
          for(int i=0 ; i<numHmmStates ; ++i)
              for(int j=0 ; j<alphabetSize ; ++j)
                  E[i][j]+=thread->E[i][j];
          for(int i=thread->first ; i<=thread->last ; ++i)
              likelihoods[i]+=thread->likelihoods[i];
          likelihoodCorrection+=thread->likelihoodCorrection;
          delete thread;
      }
      
      /********************************************************************
	Recompute the a's and e's from the A's and E's (i.e., probabilities 
	from counts)
      ********************************************************************/
      for(int k=1 ; k<numHmmStates ; ++k)
	{
	  sum=0;
	  for(Symbol s=0 ; s<alphabetSize ; ++s) sum+=E[k][s];
	  if(sum==0) cout<<String("WARNING: HMM state ")+k+
              " cannot emit any symbol; possibly unreachable."<<endl;
	  else for(Symbol s=0 ; s<alphabetSize ; ++s)
	    {
	      double newValue=E[k][s]/sum;
	      hmm.setEmissionProb(k,s,newValue);
	      if(newValue==1)
		cerr << "WARNING: hmm.getEmitProb(k,s)=1 for k=" 
		     << k << ", s=" << alphabet.lookup(s) << endl; 
	    }
	}
      for(int k=0 ; k<numHmmStates ; ++k)
	{
	  sum=0;
	  for(int ll=0 ; ll<numHmmStates ; ++ll) sum+=A[k][ll];
	  if(sum==0) for(int ll=0 ; ll<numHmmStates ; ++ll) 
	    hmm.setTransitionProb(k,ll,0);
	  else for(int l=0 ; l<numHmmStates ; ++l)
	    {
	      double p=A[k][l]/sum;
	      hmm.setTransitionProb(k,l,p);
	    }
	}

      /*********************************************************************
	Compute the log likelihood.  This is only necessary during 
        development and debugging, to ensure that recent changes to the 
        code have not changed the invariant that likelihood should increase
        monotonically (possibly except for tiny fluctuations due to 
        rounding in a digital computer):
      *********************************************************************/
      double logLikelihood=likelihoodCorrection;
      for(int j=0 ; j<numTrain ; ++j)
	if(trainingSet[j]->getLength()>0)
	  logLikelihood+=log(likelihoods[j]);
      deltaLikelihood=logLikelihood-oldLogLikelihood;
      oldLogLikelihood=logLikelihood;
      cerr << "log(likelihood)=" << logLikelihood;
      if(i>0)
	{
	  cerr << " " << progress.getProgress(i) << endl; 
	  if(deltaLikelihood<0)
	    cerr << "^--WARNING!  LIKELIHOOD DECREASED--^" << endl;
	}
      else cerr<<endl;
    }

  /*************************************************************************
   Allow start-state-self-transitions only if empty strings are to be 
   permitted
   ************************************************************************/
  if(allowEmptyStrings)
    {
      double p=double(numEmptyStrings)/double(numTrain);
      hmm.setTransitionProb(0,0,p);
      double q=1-p;
      for(int state=1 ; state<numHmmStates ; ++state)
	{
	  double rescaled=hmm.getTransitionProb(0,state)*q;
	  hmm.setTransitionProb(0,state,rescaled);
	}
    }
}



/*
  Update the expected emission counts for the HMM for a training sequence.
 */
void BWThread::reviseExpectedEmitCounts(int sequenceLength,
                                        ForwardAlgorithm &f,
                                        BackwardAlgorithm &b,
                                        Sequence &sequence)
{
  double seqProb=f.getScaledP();
  for(int k=1 ; k<states ; ++k)
    for(Symbol s=0 ; s<alphabetSize ; ++s)
      {
	double sum=0.0;
	for(int i=1 ; i<=sequenceLength ; ++i)
	  if(s==sequence[i-1])
	    {
	      sum+=f.getScalingFactors()[i] * scalingRatios[i]
		* f(k,i) * b(k,i) / seqProb;
	      if(isinf(sum) || isnan(sum)) throw String("sum=")+sum+
		  " in ParallelBaumWelch::reviseExpectedEmitCounts()";
	    }
	E[k][s]+=sum;
      }
}



/*
  Update the expected transition counts for the HMM for a training sequence.
 */
void BWThread::reviseExpectedTransCounts(int sequenceLength,
                                         ForwardAlgorithm &f,
                                         BackwardAlgorithm &b,
                                         Sequence &sequence)
{
  double seqProb=f.getScaledP();
  for(int k=0 ; k<states ; ++k)
    for(int l=1 ; l<states ; ++l)
      {
	double sum=0.0;
	for(int i=0 ; i<sequenceLength ; ++i)
	  {
	    sum+=f(k,i) * hmm.getTransitionProb(k,l) * 
	      hmm.getEmissionProb(l,sequence[i]) * b(l,i+1) / seqProb * 
	      scalingRatios[i+1];
	    if(isinf(sum) || isnan(sum)) throw String("sum=")+sum+
		   " in ParallelBaumWelch::reviseExpectedTransCounts()";
	  }
	A[k][l]+=sum;
      }
  for(int k=1 ; k<states ; ++k)
    A[k][0]+=f(k,sequenceLength) * hmm.getTransitionProb(k,0) / seqProb *
      scalingRatios[sequenceLength] * f.getScalingFactors()[sequenceLength];
  A[0][0]=0;
}



