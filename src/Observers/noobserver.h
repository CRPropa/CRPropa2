/**
   @file   noobserver.h
   @author Eric Armengaud, armengau@in2p3.fr
   @brief  Class describing no observer (in full trajectories mode)
*/

#ifndef _NOOBSERVER_H_
#define _NOOBSERVER_H_

#include "observers.h"

/**
   @class TNoObserver
   @brief Used when the full trajectories are recorded.
*/

class TNoObserver : public TObservers {
 public:
	TNoObserver() {SetType(OBSERVER_NO);}
};

#endif
