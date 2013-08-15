/**
   @file    list1dphotons.h
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Class describing the list of EM showers in a 1d environment
*/

#ifndef _LIST1DPHOTON_H_
#define _LIST1DPHOTON_H_

#include "photon.h"
#include <vector>

using namespace std;

#define USE_GRID_POSITIONS_1DPHOTONS 1 
/**< If this flag is set to 1, a predefined list of positions is used as a grid of electromagnetic shower initial positions. */
#define DEFAULT_STEP_1DPHOTONSHOWERS_MPC 5. 
/**< Stepsize of the grid of electromagnetic showers, if the flag USE_GRID_POSITIONS_1DPHOTONS is set to 0. It determines the number of showers that will be propagated. */

/**
   @class TList1DPhotons
   @brief List of electromagnetic cascades in the 1-dimensionnal case. 

   A list of empty cascades, located on fixed positions, is built at the beginning of the simulation. These cascades are then filled during the nucleon propagation, and propagated in the end.

*/

class TList1DPhotons {
 public :
  TList1DPhotons(TUniverse*, char* type="1D_ALL") ;
  ~TList1DPhotons() ;

  TPhoton* ListPhotons (int i) const { return _fListPhotons.at(i); } ;
  /**< Vector of cascades. */
  double ListInitPositions(int i) const { return _fListInitPositions.at(i); } ;
  /**< Vector of initial cascade positions. */
  int NbPhotons() const { return _fNbPhotons; } ;
  /**< Number of cascades. */

  void AddPhotons(double, double, double, int) ;
  /**< Fills a cascade with the secondaries from pair production on the CMB. The Mass is now taken into account. */
  void AddPairProdPhotons(double, double, double, const TUniverse*, double) ;
  /**< Same as above : pair production on the CMB, but with Kelner parametrization */
  void AddPhotons(double, double, PARTICLE) ;
  /**< Fills a cascade with a single particle at a given distance. */
  void AddPhoton(const TUniverse*) ;
  /**< Direct photon injection. */
  void AddPhotons(double, double, double, double, double) ;
  /**< Direct photon injection. */

  void Propagate(const TUniverse*, TBasicParam*) ;
  /**< Propagate each of the individual cascades. */

 protected :
  vector<TPhoton*> _fListPhotons ;
  vector<double> _fListInitPositions ;
  int _fNbPhotons ;
  double _fStepShowers ;

};

#endif
