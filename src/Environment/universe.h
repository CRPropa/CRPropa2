/**
   @file   universe.h
   @author Tristan Beau, beau@in2p3.fr
   @brief  Class describing the environment of the simulation : sources, targets, magnetic fields, interactions.
*/

#ifndef _UNIVERSE_H_
#define _UNIVERSE_H_

#include "sources.h"
#include "discretesources.h"
#include "continuoussources.h"

#include "magfield.h"
#include "nullmagfield.h"
#include "field1d.h"
#include "uniformmagfield.h"
#include "lssfield.h"
#include "kolmogoroffmagfield.h"

#include "observers.h"
#include "noobserver.h"
#include "pointobserver.h"
#include "smallsphereobservers.h"
#include "largesphereobservers.h"

#include "interactiondata.h"
#include "basicpinteractions.h"
#include "nullinteractions.h"
#include "sophiainteractions.h"
#include "photoninteractions.h"

#include "variableinfrared.h"
#include "uniforminfrared.h"

#include "gas.h"
#include "gridgas.h"
#include "lssgas.h"
#include "clustergas.h"

#include <fstream>
#include <iostream>
#include <vector>

#include "units.h"
#include "typedclass.h"

using namespace std;

class TIRBzEvolutionModel;

/**
   @class TUniverse
   @brief Virtual class describing all the simulation environment.

   A TUniverse always contains pointers to fundamental objects such as TSources, TMagField, TObservers,
   TInteractionData, TInfrared, and a radio background flag.
   The Bflag() function is TRUE if we are in 3D and there are magnetic fields.

*/

class TUniverse : public TTypedClass {
 public:

  virtual ~TUniverse() {}
  
  TSources* Sources() const { return _fpSources; }
  TMagField* MagField() const { return _fpMagField; }
  TObservers* Observers() const { return _fpObservers; }
  TInteractionData* InteractionData() const { return _fpInteractionData; }
  TIRBzEvolutionModel* IRBzEvolutionModel() const { return _fModel; }
  IR* Infrared() const { return _fpInfrared; }
  TGas* Gas()     const { return _fpGas; }

  bool Bflag() const { return _fBflag; }
  /**< Returns 1 if we are in 3D and there are magnetic fields (i.e. deflections must happen) */

  virtual double Xmin() const { return 0; }
  virtual double Xmax() const { return 0; }
  virtual double Ymin() const { return 0; }
  virtual double Ymax() const { return 0; }
  virtual double Zmin() const { return 0; }
  virtual double Zmax() const { return 0; }

  virtual string Integrator() const { return " "; }
  virtual double IntegratorEps() const { return 0; }
  virtual double IntegratorMinTimeStep() const { return (1.e-6 * Mpc * inv_c_light) ; }

  virtual string RadioBackground() const { return _fRadioBackground; }
  /**< Returns a flag for radio background (High, Med or Obs) useful for EM shower development */

  virtual double DistanceArray(int) const { return 0; }
    virtual std::vector<double> DistanceArray() { return std::vector<double>(); }
  virtual double RedshiftArray(int) const { return 0; }
  virtual double OmegaM() const { return 0; }
  virtual double OmegaLambda() const { return 0; }
  virtual double H0() const { return 0; }

protected:
	  
  bool _fBflag ; /**< Flag for deflections */
  
  string _fRadioBackground ; /**< Flag for the radio background */ 
  
  TSources *_fpSources ; /**< Pointer to the source(s) */
  TMagField *_fpMagField ; /**< Pointer to the magnetic field */
  TObservers *_fpObservers ; /**< Pointer to the observer(s) */
  TInteractionData *_fpInteractionData ; /**< Pointer to the interaction objects */
  IR *_fpInfrared ; /**< Pointer to the infrared background */
  TGas *_fpGas ; /**< Pointer to the gas distribution */
  TIRBzEvolutionModel* _fModel;
};
#endif
