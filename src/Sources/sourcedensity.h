/**
   @file    sourcedensity.h
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Class describing source density for discrete/continuous sources
*/

#ifndef _SOURCEDENSITY_H_
#define _SOURCEDENSITY_H_

#include "vector3d.h"
#include <math.h>
#include <fstream>
#include <vector>
#include "units.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGeneral.h"
#include "xmlparam.h"

#include "fitsio.h"

#include "fits_err.h"

/**
   @class TSourceDensity
   @brief Class containing a 1- or 3-dimensionnal array of source density, or simple grid boundaries in the case of a uniform density.

   1) Case of a 1D "Grid" type density. An ASCII array giving the COMOVING density as a function of position must be given, and positions will be drawn from this density.
   2) Case of a 3D "Grid" type density. A 3x3 array (FITS or ASCII) must be given. To draw a position from this array, a cell (i,j,k) is choosen from this array, and the position is then drawn inside this cell according to a uniform distribution.
   3) Case of a "Uniform" type density. In that case, Xmax [and optionnally Xmin] (1D case) or X/Y/Zmin/max (3D case) must be specified and the positions are drawn uniformly in the specified space region.
   If a magnetic grid is also specified, the dimensions of the boxes defined by the magnetic grid and by the source density must coincide.

*/

class TSourceDensity : public TXmlParam {
 public:
  TSourceDensity(const char*) ;
  ~TSourceDensity() ;

  bool IsUniform() const { return _fIsUniform; }
  /**< Flag set to 0 for a grid density, and to 1 for a uniform density. */
  double* DensityArray() const {return _fpDensityArray; }
  /**< Array of dimensions (Nx*Ny*Nz) containing the grid density. NULL pointer in the case of a uniform density. */
  int Nx() const { return _fNx; }
  /**< Number of cells along the x dimension for the grid density. Set to 0 for a uniform density. */
  int Ny() const { return _fNy; }
  /**< Number of cells along the y dimension for the grid density. Set to 1 in the 1-dimensionnal case. Set to 0 for a uniform density. */
  int Nz() const { return _fNz; }
  /**< Number of cells along the z dimension for the grid density. Set to 1 in the 1-dimensionnal case. Set to 0 for a uniform density. */
  double Stepsize() const { return _fStepsize; }
  /**< Stepsize of the density grid. Must be specified (in Mpc) in the configuration file. */
  TVector3D Origin() const { return _fOrigin; }
  /**< (X,Y,Z) coordinates for the origin of the density grid. Set to 0 by default. */
  TVector3D getSourcePosition() const ;
  /**< Returns a position randomly drawn from the density grid, or uniformly drawn in the case of a "Uniform" density. */
  double getDensity(TVector3D) const ;
  /**< Interpolates the density in a given place from the grid. Uses a trilinear interpolation of the grid. Works only if there is a grid. */
  double Xmin() { return _fXmin; }
  /**< Boundary of the density box. If not specified, set to 0 in the 1-D case */
  double Xmax() { return _fXmax; }
  /**< Boundary of the density box. */
  double Ymin() { return _fYmin; }
  /**< Boundary of the density box. Set to 0 in the 1-D case */
  double Ymax() { return _fYmax; }
  /**< Boundary of the density box. Set to 0 in the 1-D case */
  double Zmin() { return _fZmin; }
  /**< Boundary of the density box. Set to 0 in the 1-D case */
  double Zmax() { return _fZmax; }
  /**< Boundary of the density box. Set to 0 in the 1-D case */

 private:
  bool _fIsUniform ;
  int _fNx, _fNy, _fNz ;
  double _fStepsize ;
  TVector3D _fOrigin ;
  double *_fpDensityArray ;
  RandGeneral *_fpRandDistriX ; // Pointer to the random distribution projected along X
  vector<RandGeneral*> _fpRandDistriXY ; // Projected along (X,Y)

  double _fXmin, _fXmax ;
  double _fYmin, _fYmax ;
  double _fZmin, _fZmax ;

};

#endif
