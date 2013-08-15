/**
   @file   pointobserver.h
   @author Eric Armengaud, armengau@in2p3.fr
   @brief  Class describing a 'punctual' observer at the origin for 1D simulations
*/

#ifndef _POINTOBSERVER_H_
#define _POINTOBSERVER_H_

#include <vector>
#include "observers.h"
#include "xmlparam.h"
#include "crp_err.h"

/**
   @class TPointObserver
   @brief Used in 1-dimensionnal simulations, where a detection takes place when a particle reaches x<0
*/

class TPointObserver : public TObservers, public TXmlParam {
 public:
  TPointObserver(const char*) ;
  ~TPointObserver() {} 

  int Nb() const { return _fNb; }
  /** Set to 1 */
  TVector3D Positions(int) const { return TVector3D(); }
  /** Set to 0 */

 private:
  int _fNb ;
};

#endif
