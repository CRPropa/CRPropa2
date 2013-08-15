/**
   @file    nullinteractions.h
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Trivial class in the absence of interactions
*/

#ifndef _NULLINTERACTIONS_H_
#define _NULLINTERACTIONS_H_

#include "interactiondata.h"
#include "units.h"

#define NOINTER_MAXSTEP_MPC 10.
/**< Bound on the propagation step, even if there are no interactions. */

/**
   @class TNullInteractions
   @brief Trivial implementation of the TInteractionData class, with all interactions switched off.

*/

class TNullInteractions : public TInteractionData {
 public:
  TNullInteractions() { 
    SetType(INTERACTION_NO); 
    _fInteractionTimeStep = NOINTER_MAXSTEP_MPC * Mpc * inv_c_light ; }
  ~TNullInteractions() {}
};

#endif
