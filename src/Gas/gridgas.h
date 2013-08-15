/**
   @file    gridgas.h
   @author  Luca Maccione, luca.maccione@desy.de
   @brief   Class describing gas density fields which are stored as a 3D array with regular mesh
*/

#ifndef _GRIDGAS_H_
#define _GRIDGAS_H_

#include "gas.h"

#include <cmath>

/**
   @class TGridGas
   @brief Generic class describing a 3D gas field on a regular grid.

   Currently, the classes TLSSGas is an implementation of this class. It contains the geometrical description of the grid, as well as the interpolation routine to estimate the gas density at any position.

*/
class TGridGas : public TGas {
 public:
  TGridGas() : TGas() { SetType(GRIDGAS) ;}
  std::vector<double> GasDensity() const { return _fpGas; }
  /**< Array of gas, of size (Nx*Ny*Nz).
     As this is a 1-dimensionnal array, a 3-dimensionnal position (i,j,k) corresponds by convention to the index I = i*NyNz + j*Nz + k */
  int Nx() const { return _fNx; }
  /**< Dimension of the gas grid along the (x) axis. */
  int Ny() const { return _fNy; }
  /**< Dimension of the gas grid along the (y) axis. */
  int Nz() const { return _fNz; }
  /**< Dimension of the gas grid along the (z) axis. */
  double Stepsize() const { return _fStepsize; }
  /**< Stepsize of the grid, specified in the configuration file in Mpc. */
  TVector3D Origin() const { return _fOrigin; }
  /**< Origin position of the grid. Set to 0 by default. */
  double Density(TVector3D);
  /**< Trilinear interpolation of the gas field at any point. */
  
 protected:
  int _fNx, _fNy, _fNz ; // size of the B array
  double _fStepsize ; // distance between 2 points in the array
  TVector3D _fOrigin ; // coordinates of the origin of the array
  std::vector<double> _fpGas ; // gas density field
};

#endif
