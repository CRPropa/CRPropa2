/**
   @file    uniformmagfield.h
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Class describing a uniform magnetic field
*/

#ifndef _UNIFORMMAGFIELD_H_
#define _UNIFORMMAGFIELD_H_

#include "magfield.h"

/**
   @class TUniformMagField
   @brief Uniform 3-dimensionnal magnetic field, with a value that is set in the configuration file.
*/

class TUniformMagField : public TMagField {
 public:
  TUniformMagField(TVector3D aField) { SetType(MAGFIELD_UNIFORM); _fFieldValue = aField; }
  ~TUniformMagField() {} 
  
  TVector3D getField(TVector3D) const { return _fFieldValue; }
  /**< Returns the field value, which is the same at any position. */
  
 private:
  TVector3D _fFieldValue ;
  
};

#endif
