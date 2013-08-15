/**
   @file    Opt.h
   @author  Nils Nierstenhoefer, nierstenhoefer@physik.uni-wuppertal.de
   @brief   Returns photon density (h_bar=c=1) 
*/

#ifndef _Opt_H_
#define _Opt_H_

#include "typedclass.h"
#include "PhotonBackground.h"

class Opt : public PhotonBackground{
 public:
  Opt();
  ~Opt();
  virtual double GetPhotonDensity(double x, 
				  double y, 
				  double z, 
				  double redshift,
				  double PhotonEnergy);
};





#endif

