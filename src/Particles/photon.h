/**
   @file    photon.h
   @author  Eric Armengaud, armengau@in2p3.fr & Joerg Kulbartz jkulbart@mail.desy.de
   @brief   Class describing (secondary) photons and electromagnetic cascade
*/

#ifndef _PHOTON_H_
#define _PHOTON_H_

#include "particle.h"
#include "prop_second.h"
#include "defqueue.h"
#include "env1d.h"

#define PHOTON_DETECTOR_CROWN 1 
/**< Flag : if 1, checks for observers inside the current simulation box, but also inside the 26 nearest boxes surrounding the particle */

#define PAIRPROD_SPECTRALINDEX -7./4.
/**< Spectral index assumed for the distribution of secondary pairs from pair production by protons on the CMB. */

#define PAIRPRODCUT_EE_EP 10
/**< For the Kelner pair production spectrum: cut the spectrum for Ee > Ep/PAIRPRODCUT_EE_EP, since the analytical formula is valid only for Ee << Ep */

using namespace std;

/**
   @class TPhoton
   @brief Implementation of secondary electromagnetic cascades (photons, electrons, positrons)

   The propagation of these cascades is achieved using the DINT code. Once an observer is found, all the relevant parameters (propagation distance, cosmological parameters, injection spectra of the electromagnetic species, magnetic field along the line of sight) are given to the main DINT driver. In the end, propagated spectra in the full MeV - kZeV energy range are recorded.

*/

class TPhoton : public TParticle {
 public:
  TPhoton(double, TUniverse*, char* type="1D_ALL") ;
  /**< Trivial constructor for the 1-dimensionnal case. Builds an empty electromagnetic cascade at a given distance from the observer. */
  TPhoton(TUniverse*) ;
  /**< Gamma sources, not implemented. */
  TPhoton(TVector3D, TVector3D, double, TVector3D, double, TUniverse*, int, string, double) ;
  /**< Secondaries from pion production. The position, momentum and charge of the secondary are given (it can be a gamma, an electron or a positron), as well as the source position and injection energy. Additionally a string for output and the time of particle creation is given. */
  TPhoton(double, double, TVector3D, TVector3D, TVector3D, double, TUniverse*, int, int, string, double) ;
  /**< Secondaries from pair production. The main parameters that are given are the energy of the primary proton and its energy loss (pair production is treated as a continuous energy loss process). The position, momentum, source position and injection energy are also specified. The e+/- spectrum is built from the energy loss and assuming a power-law distribution. For nuclei the mass of the nucleus is now taken into account. Additionally a string for output and the time of particle creation is given. */

  TPhoton(TUniverse*, double, double, double, TVector3D, TVector3D, TVector3D, double, double, int, int, string, double) ; 
  /**< Same as above : secondaries from pair production; but using Kelner parametrization, so that the energy of the primary proton and its timestep are given as arguments.  Additionally a string for output and the time of particle creation is given. */
  ~TPhoton() ;

  QUEUE<TParticle*>* Propagate(TBasicParam*) ;
  /**< Driver routine for the electromagnetic cascade detection and propagation. */
  vector<double> CheckDetection(TBasicParam*) ;
  /**< Detection by an observer. */
  void DevelopShower() ;
  /**< Main routine running the interface with DINT. The relevant parameters (magnetic field along the line of sight, IR background, cosmological parameters) are built and the routine prop_second is called to propagate the initial electromagnetic spectrum on a given distance. */
  void WriteShower(TBasicParam*) const ;
  /**< Photon spectrum output routine. */

  void AddToSpectrum(double, double, int) ;
  /**< Add a power-law spectrum due to pair production to the electromagnetic cascade. The given parameters are the energy of the primary proton and its energy loss. Additionally the Mass is taken into account. */

  void AddPairProdSpectrum(double, double, double, const TUniverse*) ;
  /**< Add particles due to pair production, according to the Kelner parametrization. The given parameters are the energy of the primary proton, the time interval during which the energy loss takes place and the redshift. */

  void AddToSpectrum(double, PARTICLE) ;
  /**< Add a single particle to the electromagnetic cascade. Used for the pion production secondaries. */
  void AddToSpectrum(double, double, double, double) ;
  /**< Add to the spectrum a power-law distribution of photons, of spectral index aAlpha, between aEmin and aEmax, with total energy aEint. Used for direct photon injection. */

  Spectrum FullSpectrum() const { return _fFullSpectrum ; } ;
  /**< The spectrum of all species (not only the photons) that is used by DINT. A "Spectrum" is a DINT structure containing an array of fixed size. */
  dCVector EnergyGrid() const { return _fEnergyGrid ; } ;
  /**< The energy grid associated to FullSpectrum. A "dCVector" is a DINT structure containing an array. */

 protected:
  vector<double> _fShowerEnergy ;
  vector<double> _fShowerSpectrum ; // Only the gamma spectrum, recorded in the end
  double _fShowerPropDistance ;
  double _fDeltaE_hadron ;
  double _fInjE_hadron ;
  int _fPairProdFlag ;
  int _fInitType; 
  TVector3D _fInitPosition ;
  TVector3D _fSourcePosition ; // source position of the proton which generated the gamma
  double _fSourceEnergy ; // energy at the injection of the proton which generated the gamma (for reweighting, in 3D)
  Spectrum _fFullSpectrum ; // all the species : spectrum to be sent to DINT
  dCVector _fEnergyGrid ; // IN DINT UNITS!
  dCVector _fEnergyWidth ;
  string _fOriginStr;
  double _fInitTime;

};

#endif
