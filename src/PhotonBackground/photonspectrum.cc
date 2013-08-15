/**
   @file    photonspectrum.cc
   @author  Luca Maccione, luca.maccione@desy.de
   @brief   Implementation of the TPhotonSpectrum class. See the .h file
*/

#include "photonspectrum.h"

using namespace std;

TPhotonSpectrum::TPhotonSpectrum(vector<double> En, vector<double> sp) {
  _fNE = En.size();
  _fEnergy = vector<double>(En);
  _fEnergyDensity = vector<double>(sp);
}

double TPhotonSpectrum::Spectrum(double en) {
  // Needs refinement. Hard code the energy binning (either linear or log) ?
  int i = 0;
  while (_fEnergy[i] < en) i++;
  int j = i-1;

  double r1 = (en-_fEnergy[j])/(_fEnergy[i]-_fEnergy[j]);
  if(en<_fEnergy[0] || en> _fEnergy[_fNE-1]){
    return 0.;
  }
  else{
    return (1.0-r1)*_fEnergyDensity[j] + r1*_fEnergyDensity[i];
  }
}

