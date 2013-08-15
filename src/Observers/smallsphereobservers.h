/**
   @file   smallsphereobservers.h
   @author Eric Armengaud, armengau@in2p3.fr
   @brief  Class describing spherical observers (detection of particle goes inside sphere)
*/

#ifndef _SMALLSPHEREOBSERVERS_H_
#define _SMALLSPHEREOBSERVERS_H_

#include "vector"
#include "observers.h"
#include "xmlparam.h"
#include "crp_err.h"

#define OBS_SECURITY_RATIO 0.05
/**<  When a particle is at a distance (1+OBS_SECURITY_RATIO)*Radius from an observer, and points toward it, it is propagated inside the sphere and recorded. */
#define OBS_TCORR_INC 0.001
/**< At the last step before recording, at a propagation time T from the sphere, it is propagated in straight line during a time (1+OBS_TCORR_INC)*T */

/**
   @class TSmallSphereObservers
   @brief Small spheres surrounding an observer. An event is recorded when a particle goes through one of these observers.

The principle of the detection algorithm is the following: when a particle points toward a sphere and is susceptible to reach it, it is propagated with steps that are small enough so that it won't "go through" the observer. When the particle is very near the sphere and points toward it, it is propagated inside it in straight line, and recorded. After that, the particle continues its way out of the sphere.

*/

class TSmallSphereObservers : public TObservers, public TXmlParam {
 public:
  TSmallSphereObservers(const char *) ;
  ~TSmallSphereObservers() {}
  
  int Nb() const { return _fNb; }
  /**< Number of small sphere observers */
  TVector3D Positions(int i) const { return _fPositions.at(i); }
  /**< Positions of the centers of the spheres. */
  double Radius() const { return _fRadius; }
  /**< Radii of the spheres. */

 private:
  int _fNb ;
  vector<TVector3D> _fPositions ;
  double _fRadius ;
};

#endif
