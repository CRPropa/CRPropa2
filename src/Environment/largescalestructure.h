/**
   @file    largescalestructure.h
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Class describing all the parameters of a large-scale 3-D environment
*/

#ifndef _LSS_H_
#define _LSS_H_

#include "universe.h"
#include "nullinteractions.h"
#include "xmlparam.h"
#include "crp_err.h"

#define MINSTEP_FACTOR 0.1
/**< Tuning of the minimal stepsize of charged particles during propagation: all propagation steps are forced to be larger than (stepsize of the magnetic grid)*MINSTEP_FACTOR. It is indeed useless to achieve much smaller steps than the magnetic grid resolution can allow. */
#define GRID_MATCH_MPC 0.01
/**< Relative accuracy with which the boundaries of the 3-dimensionnal magnetic and density grids must coincide. */

/**
   @class TLargeScaleStructure
   @brief Environment for 3D simulations on long distances

   This environment is a 3D box, whose dimensions are defined by X/Y/Zmin/max, and with periodic boundary conditions assumed for CR propagation.
   The magnetic field can be a 3D grid, the sources can be continuous or discrete.
   In an "event" mode, the Observers can be either large spheres surrounding a volume where the sources stand, or small spheres that play the role of "detectors".

*/

class TLargeScaleStructure : public TUniverse, public TXmlParam {
 public:
  TLargeScaleStructure(const char*) ;
  ~TLargeScaleStructure() ;

  string Integrator() const { return _fIntegrator; }
  double IntegratorEps() const { return _fIntegratorEps; }
  double IntegratorMinTimeStep() const { return _fIntegratorMinTimeStep; }
  double Xmin() const { return _fXmin; }
  double Xmax() const { return _fXmax; }
  double Ymin() const { return _fYmin; }
  double Ymax() const { return _fYmax; }
  double Zmin() const { return _fZmin; }
  double Zmax() const { return _fZmax; }

 private:
  string _fIntegrator ;
  double _fIntegratorEps ;
  double _fIntegratorMinTimeStep ;
  double _fXmax, _fXmin ;
  double _fYmax, _fYmin ;
  double _fZmax, _fZmin ;

};

#endif
