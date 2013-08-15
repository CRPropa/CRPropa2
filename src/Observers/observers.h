/**
   @file   observers.h
   @author Eric Armengaud, armengau@in2p3.fr
   @brief  Virtual class describing any kind of observers in the simulation
*/

#ifndef _OBSERVER_H_
#define _OBSERVER_H_

#include "typedclass.h"
#include "vector3d.h"
#include "units.h"
#include <fstream>
#include <iostream>

using namespace std;

/**
   @class TObservers
   @brief Virtual class for all observers. Currently, apart from the trivial classes TNoObserver and TPointObserver (corresponding to the full propagation mode and the 1D case), there are two implementations of this class: TSmallSphereObserver and TLargeSphereObserver.

*/

class TObservers : public TTypedClass {
 public:

  virtual ~TObservers() {}

  virtual int Nb() const { return 0; }
  virtual TVector3D Positions(int) const { return TVector3D(); }
  virtual double Radius() const { return 0; }
  virtual double Radii(int) const { return 0; }

};

#endif
