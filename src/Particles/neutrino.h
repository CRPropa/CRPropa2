
/**
   @file    neutrino.h
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Class describing (secondary) neutrinos
*/

#ifndef _NEUTRINO_H_
#define _NEUTRINO_H_

#include "defqueue.h"
#include "particle.h"
#include <string.h>

#define NEUTRINO_DETECTOR_CROWNS 1 
/**< Number of crowns surrounding the propagation box where one must check for observers for neutrinos. If set to 0, only observers inside the same box as the one where the neutrino is injected will be searched. If set to 1, 3^3-1=26 more nearby boxes will be searched; if set to 2, 124 nearby boxes are searched. */

using namespace std;

/**
   @class TNeutrino
   @brief Implementation of neutrinos, simply inherited from the TParticle class.

   After they are generated as secondaries from UHECR interactions, neutrinos are checked for detection by an observer, and if the detection is possible they are propagated in one step to the observer: straigh-line propagation, only redshift effects are taken into account. No oscillation phenomena are modellized here (they can be taken into account a posteriori, from the output files).

*/

class TNeutrino : public TParticle {
 public:
  TNeutrino(int, TVector3D, TVector3D, double, TUniverse*, int, double) ;
  /**< Constructor from the flavor, position, momentum and the injection energy, type of the primary and initial time. */
  TNeutrino(int, TVector3D, TVector3D, TVector3D, double, TUniverse*, int, double) ;
  /**< Constructor from the flavor, position and momentum of the secondary neutrino, and the injection position, energy, type of the primary and time. */
  TNeutrino(TUniverse*) ;
  /**< Constructor of neutrinos from sources (not as secondaries). NOT IMPLEMENTED. */
  ~TNeutrino() {}

  QUEUE<TParticle*>* Propagate(TBasicParam*) ;
  /**< Driver routine for neutrino propagation. */
  vector<double> CheckDetection(TBasicParam*) ;
  /**< Checks detection by the observers. There are some differences w.r.t. the TParticlePropa routine, due to the fact that there is only one propagation step. */
  void Detect(TBasicParam*, vector<double>) const ;
  /**< Prints a neutrino "event" to the output. */

 private:
  int _fInitType;
  TVector3D _fInitPosition ;
  TVector3D _fSourcePosition ; // source position of the proton which generated the nu
  int _fFlavor ; // following HepPDT convention (nu_e=12,nu_e~=-12,nu_mu=14,nu_mu~=-14)
  double _fInitialEnergy ; // will change with redshift in 1D!
  double _fSourceEnergy ; // to know injection energy of the primary (for reweighting)
  double _fInitTime;
};

#endif
