#ifndef _UNIFORMINFRARED_H_
#define _UNIFORMINFRARED_H_

#include "PhotonBackground.h"
#include "photonspectrum.h"
#include "IR.h"
#include "uniforminfrared.h"

#include "xmlparam.h"
#include "xmlextract.h"
#include "crp_err.h"

#include <vector>
#include <string>
#include <fstream>
#include <math.h>

class TUniformInfrared : public IR, public TXmlParam {
 public:
  TUniformInfrared();
  TUniformInfrared(const char*) ;
  virtual ~TUniformInfrared() { delete spectrum; } 
  
  virtual double GetPhotonDensity(double x, 
				  double y, 
				  double z, 
				  double redshift,
				  double PhotonEnergy);

  virtual double Zmax() { return _fZmax; }
  /**< Maximum redshift of the IR background. */
  //  TPhotonSpectrum Spectrum(double, double, double);
TPhotonSpectrum Spectrum() { return *spectrum;}
 private:
  TPhotonSpectrum* spectrum;
};

#endif
