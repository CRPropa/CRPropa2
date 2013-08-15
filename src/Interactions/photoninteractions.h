/**
   @file    photoninteractions.h
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Class describing interactions for electromagnetic cascades only.
*/

#ifndef _PHOTONINTERACTIONS_H_
#define _PHOTONINTERACTIONS_H_

#include "interactiondata.h"
#include "units.h"
#include "xmlparam.h"
#include "crp_err.h"

#include <vector>
#include <fstream>
#include <math.h>
#include <string.h>

/**
   @class TPhotonInteractions
   @brief Class describing interactions for the electromagnetic cascades only. This class is only useful when photons only are injected at the sources.

*/

class TPhotonInteractions : public TInteractionData, public TXmlParam {
 public:
  TPhotonInteractions(const char*) ;
  ~TPhotonInteractions() ;

  string ShowerTableDir() const { return _fShowerTableDir; }
  /**< Directory where the electromagnetic interaction tables must be read by DINT. */
  double CutcascadeFlag() const { return _fCutcascadeFlag; }
  /**< If non zero, "cuts" the e+/- cascade by the magnetic deflections in DINT: only the electrons with  r_Larmor > CutcascadeFlag() * min(t_synchrotron,t_ICS) are kept */
  double PairProdSpec(double, int) const ;
  /**< Returns the interpolated pair production spectrum for a given nucleon energy, at a given electron energy bin. The spectrum is E^2 dN/dE in eV/s, where N=N(e-)=N(e+). */
 protected:
  string _fShowerTableDir ;
  double _fCutcascadeFlag ;
  double _fPairProdTable[NBINS_PAIRPROD][NUM_MAIN_BINS];
  double _fNucleonEnergy_PP[NBINS_PAIRPROD];

};

#endif
