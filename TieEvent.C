/****************************************************************
 TieEvent.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/Exceptions.H"
#include "BOOM/HashTable.H"
#include "TieEvent.H"
using namespace std;
using namespace BOOM;


struct EntityNames {
  EntityNames() {
    hashish("means")=TIE_MEANS;
    hashish("mean")=TIE_MEANS;
    hashish("covariance_matrix")=TIE_COV_MATRIX;
    hashish("covariance_matrices")=TIE_COV_MATRIX;
    hashish("variances")=TIE_VARIANCES;
    hashish("correlation_matrix")=TIE_COR_MATRIX;
    hashish("correlation_matrices")=TIE_COR_MATRIX;
    hashish("weights")=TIE_MIXTURE_WEIGHTS;
    hashish("transitions")=TIE_TRANSITIONS;
    hashish("transition")=TIE_TRANSITIONS;
    hashish("chains")=TIE_MARKOV_CHAINS;
    hashish("chain")=TIE_MARKOV_CHAINS;
  }
  TieEntity &hashish(const char *s) {return hash[String(s)];}
  static HashMap<String,TieEntity> hash;
};
HashMap<String,TieEntity> EntityNames::hash(499);
EntityNames global;


/****************************************************************
                         TieEntity methods
 ****************************************************************/
TieEntity stringToTieEntity(const String &string)
{
  if(EntityNames::hash.isDefined(string)) return EntityNames::hash[string];
  throw string+" : unknown entity in tie profile";
}



bool canTieStates(TieEntity e)
{
  switch(e) 
    {
    case TIE_MEANS:
    case TIE_COV_MATRIX:
    case TIE_VARIANCES:
    case TIE_COR_MATRIX:
    case TIE_TRANSITIONS:
      return false;
    case TIE_MIXTURE_WEIGHTS:
    case TIE_MARKOV_CHAINS:
      return true;
    }
  INTERNAL_ERROR;
}



bool canTieComponents(TieEntity e)
{
  switch(e) 
    {
    case TIE_MEANS:
    case TIE_COV_MATRIX:
    case TIE_VARIANCES:
    case TIE_COR_MATRIX:
    case TIE_TRANSITIONS:
      return true;
    case TIE_MIXTURE_WEIGHTS:
    case TIE_MARKOV_CHAINS:
      return false;
    }
  INTERNAL_ERROR;
}



bool canFixInStates(TieEntity e)
{
  switch(e) 
    {
    case TIE_MEANS:
    case TIE_COV_MATRIX:
    case TIE_VARIANCES:
    case TIE_COR_MATRIX:
    case TIE_TRANSITIONS:
      return false;
    case TIE_MIXTURE_WEIGHTS:
    case TIE_MARKOV_CHAINS:
      return true;
    }
  INTERNAL_ERROR;
}




/****************************************************************
                         TieLevel methods
 ****************************************************************/
TieLevel stringToTieLevel(const String &s)
{
  if(s=="states") return TIE_STATES;
  if(s=="state") return TIE_STATES;
  if(s=="components") return TIE_COMPONENTS;
  throw s+" : unknown level in tie profile";
}



/****************************************************************
                         TieOp methods
 ****************************************************************/
TieOp stringToTieOp(const String &s)
{
  if(s=="tie") return TIE_PARMS;
  if(s=="fix") return FIX_PARMS;
  throw s+" : unknown operation in tie profile";
}



/****************************************************************
                         TieEvent methods
 ****************************************************************/
TieEvent::TieEvent()
  : entity(TIE_ENTITY_NONE)
{
  // ctor
}



TieEvent::TieEvent(TieOp op,TieEntity e)
  : op(op), entity(e)
{
  // ctor
}



TieOp TieEvent::getOp() const
{
  return op;
}



void TieEvent::addState(STATE s)
{
  stateList.push_back(s);
}



TieEntity TieEvent::getEntity() const
{
  return entity;
}



const BOOM::Vector<STATE> &TieEvent::getStates() const
{
  return stateList;
}



/****************************************************************
                    TieTransitionEvent methods
 ****************************************************************/
TieTransitionEvent::TieTransitionEvent(TieOp op)
  : TieEvent(op,TIE_TRANSITIONS)
{
  // ctor
}



void TieTransitionEvent::addTransition(STATE from,STATE to)
{
  addState(from);
  addState(to);
}



int TieTransitionEvent::getNumTransitions() const
{
  return getStates().size()/2;
}



void TieTransitionEvent::getTransition(int index,STATE &from,STATE &to) const
{
  const BOOM::Vector<STATE> &states=getStates();
  from=states[index*2];
  to=states[index*2+1];
}




