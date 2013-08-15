/**
   @file    gridfield.h
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Class describing magnetic fields which are stored as a 3D array with regular mesh
*/

#ifndef _GRIDFIELD_H_
#define _GRIDFIELD_H_

#include "magfield.h"
#include "fitsio.h"
#include <math.h>

/**
   @class TGridField
   @brief Generic class describing a 3D magnetic field on a regular grid.

   Currently, the classes TKolmogoroffMagField and TLSSField are implementations of this class. It contains the geometrical description of the grid, as well as the interpolation routine to estimate the magnetic field at any position.

*/
class TGridField : public TMagField {
 public:
  TGridField() { SetType(MAGFIELD_GRID) ;}
  TVector3D* B() const { return _fpB; }
  /**< Array of magnetic field (a vector), of size (Nx*Ny*Nz).
     As this is a 1-dimensionnal array, a 3-dimensionnal position (i,j,k) corresponds by convention to the index I = NyNz + j*Nz + k */
  int Nx() const { return _fNx; }
  /**< Dimension of the magnetic grid along the (x) axis. */
  int Ny() const { return _fNy; }
  /**< Dimension of the magnetic grid along the (y) axis. */
  int Nz() const { return _fNz; }
  /**< Dimension of the magnetic grid along the (z) axis. */
  double Stepsize() const { return _fStepsize; }
  /**< Stepsize of the grid, specified in the configuration file in Mpc. */
  TVector3D Origin() const { return _fOrigin; }
  /**< Origin position of the grid. Set to 0 by default. */
  TVector3D getField(TVector3D) const ;
  /**< Trilinear interpolation of the magnetic field at any point. */
  int writeMagField(string);
  /**< Writes the magnetic field to the given filename. */

 protected:
  int _fNx, _fNy, _fNz ; // size of the B array
  double _fStepsize ; // distance between 2 points in the array
  TVector3D _fOrigin ; // coordinates of the origin of the array
  TVector3D *_fpB ; // magnetic field
};

#endif
