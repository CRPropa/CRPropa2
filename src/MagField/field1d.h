/**
   @file    field1d.h
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Class describing a perpendicular magnetic field for 1D propagation of electromagnetic showers
*/

#ifndef _FIELD1D_H_
#define _FIELD1D_H_

#include "magfield.h"
#include "xmlparam.h"
#include "crp_err.h"

#include <vector>
#include <math.h>

#define CHECK_STEP 1.e-3
/**< Relative accuracy with which the stepsize of the magnetic field grid must be constant. */

/**
   @class TField1D
   @brief Magnetic field perpendicular to the line of sight in the context of 1D simulations.

   This object contains a 1D grid of field magnitudes that allows to compute the sychrotron radiation of electromagnetic secondaries at any position along the line of sight.

*/

class TField1D : public TMagField, public TXmlParam {
 public:
  TField1D(const char*) ;
  ~TField1D() {}

  double getPositions(int i) const { return _fPositions.at(i); }
  /**< Allows to parse the 1D grid of positions corresponding to the tabulated field. */
  double getFieldValues(int i) const { return _fFieldValues.at(i); }
  /**< Allows to parse the 1D grid of magnetic field */
  double Xmax() const { return _fXmax; }
  /**< Maximum distance until which the field is tabulated. */
  double Stepsize() const { return _fStepsize; }
  /**< Stepsize for the tabulated magnetic field. */

  TVector3D getField(TVector3D) const ;
  /**< Linear interpolation of the tabulated magnetic field at any position. The (x) and (z) components are null, so that B = B_y is perpendicular to the line of sight (only this component is useful to estimate the synchrotron radiation). */

 private: 
  vector<double> _fPositions ;
  vector<double> _fFieldValues ;
  double _fXmax ;
  double _fStepsize ;

};

#endif
