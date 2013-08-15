/**
   @file    lssfield.h
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Class describing magnetic fields from a LSS simulation output
*/

#ifndef _LSSFIELD_H_
#define _LSSFIELD_H_

#include "gridfield.h"
#include "xmlparam.h"
#include "crp_err.h"

#include "fitsio.h"

#include "fits_err.h"

/**
   @class TLSSField
   @brief Grid magnetic field that is read from an ASCII or a FITS file

   Magnetic field given by an external regular grid, for example derived from the LSS simulation frameworks
*/

class TLSSField : public TGridField, public TXmlParam {
 public:
  TLSSField(const char*) ;
  ~TLSSField() ;
  double Bmax() const { return _fBmax; }
  /**< Maximum value of the field on the grid. */

 private:
  double _fBmax ;
};

#endif
