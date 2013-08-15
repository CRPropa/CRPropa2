
/**
   @file   largesphereobservers.h
   @author Eric Armengaud, armengau@in2p3.fr
   @brief  Class describing spherical observers around a source
*/

#ifndef _LARGESPHEREOBSERVERS_H_
#define _LARGESPHEREOBSERVERS_H_

#include "observers.h"
#include "xmlparam.h"
#include "crp_err.h"
#include "vector3d.h"
#include <stdlib.h>
#include "defqueue.h"
#include "units.h"
#include <vector>

/**
   @class TLargeSphereObservers
   @brief Large spheres around a source.

   Particles are initially inside these large spheres, and they are recorded when they go through one of those spheres.

*/

class TLargeSphereObservers : public TObservers, public TXmlParam {
 public:
  TLargeSphereObservers(const char *) ;
  ~TLargeSphereObservers() {}

  int Nb() const { return _fNb; }
  /**< Number of spheres. */
  TVector3D Positions(int i) const { return _fPositions[i]; }
  /**< Positions of the sphere centers */
  double Radii(int i) const { return _fRadii[i]; }
  /**< Sphere radii. */

 private:
  int _fNb ;
  vector<TVector3D> _fPositions ;
  vector<double> _fRadii ;
};

#endif
