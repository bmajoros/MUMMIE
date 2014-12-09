/****************************************************************
 TieEvent.H
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_TieEvent_H
#define INCL_TieEvent_H
#include <iostream>
#include "BOOM/Vector.H"
#include "BOOM/String.H"
#include "State.H"
using namespace std;
using namespace BOOM;

/*
 Examples:
  "tie means between components"
  "tie covariance_matrix between components"
  "tie variances between components"
  "tie correlation_matrix between components"
  "tie weights between states 1,2,3"
  "tie transitions 1->2 3->4 5->6"
  "tie chains between states 5,7,12"
  "fix means"
  "fix covariance_matrix"
  "fix weights in state 8"
  "fix transition 3->6"
  "fix chain in state 7"
 */

/****************************************************************
                          enum TieEntity
 ****************************************************************/
enum TieEntity {
  TIE_ENTITY_NONE,
  TIE_MEANS,
  TIE_COV_MATRIX,
  TIE_VARIANCES,
  TIE_COR_MATRIX,
  TIE_MIXTURE_WEIGHTS,
  TIE_TRANSITIONS,
  TIE_MARKOV_CHAINS
};
TieEntity stringToTieEntity(const String &);
bool canTieStates(TieEntity);
bool canTieComponents(TieEntity);
bool canFixInStates(TieEntity);

/****************************************************************
                           enum TieLevel
 ****************************************************************/
enum TieLevel {
  TIE_LEVEL_NONE,
  TIE_STATES,
  TIE_COMPONENTS
};
TieLevel stringToTieLevel(const String &);

/****************************************************************
                           enum TieOp
 ****************************************************************/
enum TieOp {
  TIE_PARMS,
  FIX_PARMS
};
TieOp stringToTieOp(const String &);

/****************************************************************
                          class TieEvent
 ****************************************************************/
class TieEvent {
public:
  TieEvent();
  TieEvent(TieOp,TieEntity);
  void addState(STATE);
  TieOp getOp() const;
  TieEntity getEntity() const;
  const BOOM::Vector<STATE> &getStates() const;
protected:
  TieOp op;
  TieEntity entity;
  BOOM::Vector<STATE> stateList;
};

/****************************************************************
                      class TieTransitionEvent
 ****************************************************************/
class TieTransitionEvent : public TieEvent {
public:
  TieTransitionEvent(TieOp);
  void addTransition(STATE from,STATE to);
  int getNumTransitions() const;
  void getTransition(int index,STATE &from,STATE &to) const;
};

#endif

