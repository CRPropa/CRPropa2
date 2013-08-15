/**
   @file    env1d.h
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Class describing all the environment for a 1-D simulation: no magnetic field, observer at the origin.
*/

#ifndef _ENV1D_H_
#define _ENV1D_H_

// Default cosmological parameters
#define DEFAULT_OMEGA_M 0.3 /**< Default value for _fOmegaM */
#define DEFAULT_OMEGA_LAMBDA 0.7 /**< Default value for _fOmegaLambda */
#define DEFAULT_H_0_KM_S_MPC 71. /**< Default value for _fH0 */

#include "universe.h"
#include "xmlparam.h"

/**
   @class TEnv1D
   @brief Environment in 1D simulations

   The environment in 1D is simply defined by cosmological parameters that are taken into account for all redshift effects as well as the evolution of photon backgrounds, and a maximum distance from the observer. 

   The sources can be discrete or continuous. In the case where only events are recorded, the observer is a simple point at the origin, so the "detection" algorithms are trivial.
   An inhomogeneous magnetic field can be implemented : its influence is taken into account for the EM cascades.

*/

class TEnv1D : public TUniverse, public TXmlParam {
 public:
  TEnv1D(const char*) ;
  /**< Initializes all the environment parameters from the XML file */
  ~TEnv1D() ;
  
  double Xmax() const { return _fXmax; }
  /**< Maximum distance allowed from the observer : injection won't happen further */
  double DistanceArray(int i) const { return _fDistanceArray.at(i); }
    std::vector<double> DistanceArray() { return _fDistanceArray; }
  /**< Array of distances, built once for all at the beginning and that allows redshift-distance conversion
   The "distance" here is the angular distance, i.e. comoving distance/(1+z) */
  double RedshiftArray(int i) const { return _fRedshiftArray.at(i); }
  /**< Array of redshifts, built once for all at the beginning and that allows redshift-distance conversion */
  double OmegaM() const { return _fOmegaM; } /**<  cosmological matter density */
  double OmegaLambda() const { return _fOmegaLambda; } /**< cosmological dark energy */
  double H0() const { return _fH0; } /**< Hubble's constant at z=0 */

 private:
  double _fXmax ;
  double _fOmegaM ;
  double _fOmegaLambda ;
  double _fH0 ;
  vector<double> _fDistanceArray ; /**< Distance = angular distance = comoving distance/(1+z) */
  vector<double> _fRedshiftArray ;

};

#endif
