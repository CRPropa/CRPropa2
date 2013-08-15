/**
   @file    lssgas.h
   @author  Luca Maccione, luca.maccione@desy.de
   @brief   Class describing gas density field from a LSS simulation output
*/

#ifndef _LSSGAS_H_
#define _LSSGAS_H_

#include "gas.h"
#include "gridgas.h"
#include "xmlparam.h"
#include "crp_err.h"

#include "fitsio.h"

#include "fits_err.h"

/**
   @class TLSSGas
   @brief Grid gas field that is read from a FITS file

   Gas density given by an external regular grid, for example derived from the LSS simulation frameworks
*/

class TLSSGas : public TGridGas, public TXmlParam {
 public:
  TLSSGas(const char*) ;
  ~TLSSGas() { }
  double Gasmax() const { return _fGasmax; }
  /**< Maximum value of the field on the grid. */

 private:
  double _fGasmax ;
};

#endif
