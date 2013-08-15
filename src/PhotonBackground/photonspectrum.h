/**
   @file    photonspectrum.h
   @author  Luca Maccione, luca.maccione@desy.de
   @brief   Class describing a photon spectrum
*/

#ifndef _PHOTONSPECTRUM_H_
#define _PHOTONSPECTRUM_H_

#include "crp_err.h"

#include <vector>
#include <string>
#include <fstream>
#include <math.h>

/**
   @class TPhotonSpectrum
   @brief Photon Spectrum

*/

using namespace std;

class TPhotonSpectrum {
 public:
  TPhotonSpectrum() {}
  TPhotonSpectrum(vector<double> En, vector<double> sp) ;
  TPhotonSpectrum(const TPhotonSpectrum& irsp) :
    _fNE(irsp._fNE),
    _fEnergy(irsp._fEnergy),
    _fEnergyDensity(irsp._fEnergyDensity)
      {}
  ~TPhotonSpectrum() {} 
  
  double Epsmin() const { return _fEnergy.front(); }
  double Epsmax() const { return _fEnergy.back(); }
  double Spectrum(double /* eV */); /* eV/cm^3 */

  vector<double> GetE() const { return _fEnergy; }
  vector<double> Spectrum() const { return _fEnergyDensity; }

  int NE() const { return _fNE; }

 private:
  int _fNE;
  vector<double> _fEnergy;
  vector<double> _fEnergyDensity;
};
#endif
