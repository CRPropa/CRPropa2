/**
   @file   sources.h
   @author Eric Armengaud, armengau@in2p3.fr
   (changed by Nils Nierstenhoefer in September 2008 to include fixed mass & charge numbers.)
   @brief  Virtual class describing any kind of sources in the simulation
*/

#ifndef _SOURCES_H_
#define _SOURCES_H_

#include "vector3d.h"
#include <iostream>
#include "units.h"
#include "typedclass.h"
#include "sourcedensity.h"
#include <vector>
using namespace std;

/**
   @class TSources
   @brief Virtual class describing the sources of charged particles.

   All the source configurations are described by their spectra and their positions. In particular, the TSources class contains a pointer to a TSourceDensity object, which can contain a 1-D or 3-D grid of density in the simulation box.
   There are two implementations of the TSources class, depending on whether the sources are continuous (TContinuousSources) or discrete (TDiscreteSources).
*/

class TSources : public TTypedClass {
 public:
  TSources() {};
  virtual ~TSources() {}
  
  vector<int> GetInitialMassNumber() const { return _fInitialMassNumber; }       //Changed by NN
  /**< Returns the mass of the nuclei which are emitted from the source. */
  vector<int> GetInitialChargeNumber() const { return _fInitialChargeNumber; }   //Changed by NN
  /**< Returns the mass of the nuclei which are emitted from the source. */
  vector<double> GetNucleiAbundance() const {return _fAbundance; }
  int SpectrumFlag() const { return _fSpectrumFlag; }
  /**< Returns 1 for a monochromatic source, and 2 for a power-law spectrum */
  double Alpha() const { return _fAlpha; }
  /**< Spectral index in the case of a power-law spectrum (otherwise set to 0) */
  
  virtual double SigAlpha() const { return 0; }
  double Ecut() const { return _fEcut; }
  /**< Case of a power-law spectrum: Cut-off energy of the spectrum. Case of a monochromatic spectrum: injection energy. */ 
  double Emin() const { return _fEmin; }
  /**< Case of a power-law spectrum: Minimum injection energy. It is set by the minimum propagation energy of charged particles, and must be larger than Ecut in any case. */
  bool IsPhotonSources() const { return _fIsPhotonSources; }
  /**< Case of direct photon injection from the sources. If this flag is set, we must be in one-dimensional mode, and the main loop is not executed. */

  virtual unsigned int Nb() const { return 0; }
  virtual TVector3D Positions(int) const { return TVector3D(); }
  virtual double EcutList(int) const { return 0; }
  virtual double EminList(int) const { return 0; }
  virtual double AlphaList(int) const { return 0; }
  virtual int RigidityFlag(void) const {return _fRigidityFlag;}
  virtual TVector3D GetPosition() const { return TVector3D(); }
  TSourceDensity* Density() const { return _fpSourceDensity; }
  /**< Pointer to the source density, if any. */

 protected:
  std::vector<int> _fInitialMassNumber ;              //Changed by NN
  std::vector<int> _fInitialChargeNumber ;            //Changed by NN
 std::vector<double> _fAbundance ;            //Changed by NN
  int _fNSpecies;
  int _fSpectrumFlag ;
  double _fAlpha ;
  //  double _fSigAlpha ; in fact only for discrete sources
  double _fEcut ;
  double _fEmin; 
  TSourceDensity *_fpSourceDensity ;
  bool _fIsPhotonSources ;
  bool _fRigidityFlag;
};

#endif
